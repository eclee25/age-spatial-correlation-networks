/*----SQL TEMPLATE----
Author: Elizabeth Lee
Date: 9/16/14	
Function: Export list of unique zip3s in SDI database

Command Line: mysql < unique_zip3.sql | sed 's/\t/,/g' > /home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/SQL_export/unique_zip3.csv
Data: flu table: SDI
*/

SELECT distinct PATIENT_ZIP3 from flu
WHERE PATIENT_ZIP3 <> 'TOT'
;






