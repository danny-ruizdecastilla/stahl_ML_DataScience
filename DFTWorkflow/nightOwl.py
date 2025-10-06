import subprocess
import time

def submit_job(script_path):
    # Submit the job
    result = subprocess.run(["sbatch", script_path], capture_output=True, text=True)
    print(result.stdout)
    
    # Extract job ID
    job_id = result.stdout.strip().split()[-1]
    return job_id

def check_job_status(job_id):
    result = subprocess.run(["squeue", "--job", job_id], capture_output=True, text=True)
    return result.stdout

def cancel_job(job_id):
    subprocess.run(["scancel", job_id])

# Example usage
job_id = submit_job("my_script.slurm")
print(f"Submitted job ID: {job_id}")

# Wait and check
time.sleep(5)
print(check_job_status(job_id))

