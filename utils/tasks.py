from celery import Celery, current_task
import subprocess
from concurrent.futures import ThreadPoolExecutor
import pysam

celery_app = Celery('tasks', broker='redis://localhost:6379/0')

@celery_app.task(bind=True, name='tasks.run_mapdamage_task')
def run_mapdamage_task(self, bam_path, ref_path, output_dir):
    try:
        # Run mapDamage2
        cmd = [
            'mapDamage',
            '-i', bam_path,
            '-r', ref_path,
            '-d', output_dir
        ]
        subprocess.run(cmd, check=True)
    except Exception as e:
        print(f"Error running mapDamage2: {e}")
        self.update_state(state='FAILURE', meta={'exc': str(e)})
        raise e


@celery_app.task(name='tasks.run_centrifuge_task')
def run_centrifuge_task(file_path, output_path, report_path, db_path):
    centrifuge_cmd = [
        "centrifuge",
        "-x", db_path,
        "-U", file_path,
        "-S", output_path,
        "--report-file", report_path
    ]
    try:
        subprocess.run(centrifuge_cmd, check=True)
    except subprocess.CalledProcessError as e:
        raise e

@celery_app.task(bind=True, name='tasks.load_bam_file')
def load_bam_file(self, bam_path):
    try:
        total_reads = 0
        processed_reads = 0

        # Get the total read count (used to track progress)
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            total_reads = sum(1 for _ in bam_file)
        
        def process_chunk(reads):
            """ Dummy function to process chunk of reads """
            nonlocal processed_reads
            for read in reads:
                # Process each read here as needed
                processed_reads += 1
                progress = (processed_reads / total_reads) * 100
                current_task.update_state(state='PROGRESS', meta={'progress': progress})

        # Use ThreadPoolExecutor to speed up processing
        with pysam.AlignmentFile(bam_path, "rb") as bam_file, ThreadPoolExecutor(max_workers=4) as executor:
            reads = []
            chunk_size = 1000  # Number of reads to process at once

            for read in bam_file:
                reads.append(read)
                if len(reads) >= chunk_size:
                    executor.submit(process_chunk, reads)
                    reads = []

            # Process remaining reads
            if reads:
                executor.submit(process_chunk, reads)

    except Exception as e:
        self.update_state(state='FAILURE', meta={'exc': str(e)})
        raise e
