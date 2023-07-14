import os
import requests
import time
import pandas as pd
from tqdm import tqdm
import zipfile


df = pd.read_csv('./pdb2sp.csv')

sleep_retry = 15
num_retry = 5
save_dir = "./pdb_files"

def zip_directory(directory_path, zip_file_path):
    with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(directory_path):
            for file in files:
                file_path = os.path.join(root, file)
                zipf.write(file_path, os.path.relpath(file_path, directory_path))

import requests
import time

def download_pdb_file(pdb_id, save_path, max_retries=3, sleep_time_in_sec=15):
    # Define the URL to download the PDB file
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    # Initialize a counter for retry attempts
    retry_count = 0

    # Retry loop
    while retry_count < max_retries:
        try:
            # Send a GET request to the URL and save the file
            response = requests.get(url)
            if response.status_code == 200:
                with open(save_path, "wb") as file:
                    file.write(response.content)
                return True
            elif response.status_code == 404:
                print(f"PDB file _ {pdb_id}.pdb not found. Exiting.")
                return False
            else:
                print(f"Failed to download PDB file. (Attempt {retry_count + 1})")
        except requests.exceptions.RequestException as e:
            print(f"Error occurred: {e}")

        retry_count += 1
        time.sleep(sleep_time_in_sec)  # Wait for sleep_time_in_sec seconds before retrying

    print(f"Exceeded maximum retry attempts. Failed to download PDB file.")
    return False


# Check if the directory already exists
if not os.path.exists(save_dir):
    # Create the directory
    os.makedirs(save_dir)
    print("Save Directory created successfully.")
else:
    print("Save Directory already exists.")


for index, row in tqdm(df.iterrows() , total=len(df)):
    save_file_path = os.path.join(save_dir , row['PDB']+"_"+ row['SP'] +'.pdb')
    pdb_id = row['PDB']
    if os.path.exists(save_file_path):
        print(f"File {save_file_path} exists.")
    else:
        download_pdb_file(pdb_id , save_file_path , num_retry , sleep_retry)
        
zip_directory(save_dir, './pdb_files.zip')