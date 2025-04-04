import paramiko
import os
import getpass

# Define server details
remote_server = 'server.example.com'  # Replace with your server address
username = 'lihao'
password = getpass.getpass('Enter password: ')

# Get the local directory and remote directory from the user
# local_directory = input('Enter the local directory to transfer: ')
# remote_directory = input('Enter the remote directory: ')
local_directory = os.path.dirname(os.path.abspath(__file__))
remote_directory = os.path.join('/scratch/lihao/', input('Enter the new remote directory name: '))

# Create SSH client
ssh = paramiko.SSHClient()
ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(remote_server, username=username, password=password)

# Create SFTP client
sftp = ssh.open_sftp()

# Function to recursively upload a directory
def upload_directory(local_dir, remote_dir):
    if not os.path.isdir(local_dir):
        raise ValueError(f"{local_dir} is not a directory")
    
    for root, dirs, files in os.walk(local_dir):
        for dirname in dirs:
            local_path = os.path.join(root, dirname)
            relative_path = os.path.relpath(local_path, local_dir)
            remote_path = os.path.join(remote_dir, relative_path)
            try:
                sftp.mkdir(remote_path)
            except IOError:
                pass  # Directory already exists
        
        for filename in files:
            local_path = os.path.join(root, filename)
            relative_path = os.path.relpath(local_path, local_dir)
            remote_path = os.path.join(remote_dir, relative_path)
            sftp.put(local_path, remote_path)

# Ensure the remote directory exists
try:
    sftp.mkdir(remote_directory)
except IOError:
    pass  # Directory already exists

# Upload the directory
print(f'Send {local_directory}/* -> {remote_directory}')
upload_directory(local_directory, remote_directory)
# Close the SFTP client
sftp.close()

# Submit the job using SLURM
slurm_script = os.path.join(remote_directory, 'sba.sh')
exe = os.path.join(remote_directory, 'abr')
ssh.exec_command(f'chmod +x {slurm_script}')
ssh.exec_command(f'chmod +x {exe}')
stdin, stdout, stderr = ssh.exec_command(f'cd {remote_directory} && sbatch {slurm_script}')
print(stdout.read().decode())
print(stderr.read().decode())

# Close the SSH connection
ssh.close()
