import os
import logging
import logging.config

log_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logging.conf')
logging.config.fileConfig(log_file_path , disable_existing_loggers=False)
# logging.config.fileConfig('logging.conf', disable_existing_loggers=False)
logger = logging.getLogger(__name__)

def fastq_autoscan(read_file):
    """Scanning readfile to determine filetype"""
    logger.debug("Scanning readfile to determine filetype")
    with open(read_file, "r", encoding="utf-8") as rfl:
        line_ct = 0
        fastq = False
        while line_ct < 100:
            for line in rfl:
                line_ct += 1
                if line.startswith("@"):
                    fastq = True
                    break
    return fastq