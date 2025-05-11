# -*- coding: utf-8 -*-
# Author: Xuanming Liu
# Date: 2025.04.13
# Annotation: define a logger for logging. The logger can record the output of a function into a log file and print it to the console in the multiprocessing.

# =============================================================================
# Packages
# =============================================================================
import os
import logging

# ==============================================================================
# Define the logger for logging
# ==============================================================================
def define_logger(log_path, task_name):
    logger = logging.getLogger(task_name)  # 使用唯一的日志记录器名称
    if not logger.handlers:  # 避免重复添加处理器
        log_filename = os.path.join(
            log_path,
            f"{task_name}.log"
        )
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
        # File handler
        file_handler = logging.FileHandler(log_filename)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        # Stream handler
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)
    return logger

if __name__ == "__main__":
    log_path = "Results/Table1_Preprocess/cv_results"
    task_name = "test"
    logger = define_logger(log_path, task_name)
    logger.info("This is a test log message.")