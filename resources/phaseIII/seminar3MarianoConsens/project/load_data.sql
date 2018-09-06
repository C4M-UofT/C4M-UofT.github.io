\COPY admissions FROM 'admissions.csv' DELIMITER ',' CSV header; 
\COPY patients FROM 'patients.csv' DELIMITER ',' CSV header; 
\COPY icustays FROM 'icustays.csv' DELIMITER ',' CSV header; 
\COPY diagnoses_icd FROM 'diagnoses_icd.csv' DELIMITER ',' CSV header; 
\COPY d_icd_diagnoses FROM 'D_ICD_DIAGNOSES.csv' DELIMITER ',' CSV header; 
