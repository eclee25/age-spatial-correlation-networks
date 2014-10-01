/*----SQL TEMPLATE----
Author: Elizabeth Lee
Date: 10/1/14
Function: Export population size for small age groups for all service places by zip3 

(Age-specific population size and any visits should be exported separate from ILI in order to deal with lack of reporting if zip3-smallAgeGroup-week combination reported no ILI. See Evernote: Age-Spatial Corr Networks, October 8, 2014 note)

Command Line: mysql < allyears_zip3_popAge.sql | sed 's/\t/,/g' > ../../SQL_export/allyears_zip3_popAge.csv
Data: flu table: SDI
*/

/* all ages population */
SELECT YEAR(flu.WEEK), flu.PATIENT_ZIP3, flu.AGEGROUP, flu.POPSTAT from flu 
WHERE flu.SERVICE_PLACE = "TOTAL" and flu.PATIENT_ZIP3 <> "TOT" and flu.POPSTAT > 0 
GROUP BY YEAR(flu.WEEK), flu.PATIENT_ZIP3, flu.AGEGROUP
;
