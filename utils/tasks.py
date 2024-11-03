from celery import Celery
import subprocess

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
