def upload_file(file, creds, overwrite=True):
    from google.cloud import storage
    from google.cloud import bigquery
    
    # Add extension to filename
    path = file + ".csv"

    # Upload to google cloud storage
    client = storage.Client.from_service_account_json(creds)
    bucket = client.get_bucket('vie_ia2030')
    blob = bucket.blob(path)
    blob.chunk_size = 5 * 1024 * 1024 # Set 5 MB blob size to prevent timeout
    blob.upload_from_filename(path)

    # Migrate to bigquery
    client = bigquery.Client.from_service_account_json(creds)
    uri = "gs://vie_ia2030/" + path
    table_id = "vaccine-impact.data." + file
    job_config = bigquery.LoadJobConfig(
        autodetect=True,
        skip_leading_rows=1,
        # The source format defaults to CSV, so the line below is optional.
        source_format=bigquery.SourceFormat.CSV)
    if overwrite:
        job_config.write_disposition = bigquery.WriteDisposition.WRITE_TRUNCATE

    load_job = client.load_table_from_uri(uri, table_id, job_config=job_config)
    load_job.result()  # Waits for the job to complete.

    destination_table = client.get_table(table_id)  # Make an API request.
    print("Loaded {} rows".format(destination_table.num_rows))

if __name__ == "__main__":
    import sys
    upload_file(sys.argv[1], sys.argv[2])

