/*----SQL TEMPLATE----
Author: Elizabeth Lee
Date: 9/11/14
Function: Export flu season week data to calculate attack rates by season for all age groups and service places
Flu season = weeks 40 to 20

export one season at a time

Command Line: mysql < fluseasonweeks_zip3_AR.sql | sed 's/\t/,/g' > outputfile.csv
Data: flu table: SDI
*/

SELECT season_cdcwk.SMALL_SEAS_NUM, flu.WEEK, season_cdcwk.WK_NUM, flu.PATIENT_ZIP3, sum(flu.ILI_m), sum(flu.ANY_DIAG_VISIT_CT), flu.POPSTAT from flu 
LEFT JOIN season_cdcwk USING (WEEK) 
WHERE flu.AGEGROUP = "TOTAL" and flu.SERVICE_PLACE = "TOTAL" and flu.PATIENT_ZIP3 <> "TOT" and flu.POPSTAT > 0 and season_cdcwk.SMALL_SEAS_NUM = 2
GROUP BY season_cdcwk.SMALL_SEAS_NUM, flu.WEEK, flu.PATIENT_ZIP3
;





