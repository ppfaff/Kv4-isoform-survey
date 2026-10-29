@echo off
cd /d %~dp0
set PYTHONDONTWRITEBYTECODE=1
D:\anaconda3\envs\hh\python.exe -u add_kv4_survey_document.py > add_kv4_survey_document_log.txt 2>&1
type add_kv4_survey_document_log.txt
pause
