/*----SQL TEMPLATE----
Author: Elizabeth Lee
Date: 10/1/14
Function: Export flu season week data to calculate incidence rates by week for small age groups for all service places
Flu season = weeks 40 to 20.

export one season at a time

(Age-specific population size and any visits should be exported separately in order to deal with lack of reporting if zip3-smallAgeGroup-week combination reported no ILI. See Evernote: Age-Spatial Corr Networks, October 8, 2014 note)

Command Line: mysql < fluseasonweeks_zip3_incidAge.sql | sed 's/\t/,/g' > /home/elee/Dropbox/Elizabeth_Bansal_Lab/SDI_Data/age_spatial_correlation_networks/SQL_export/fluseasonweeks_zip3_incidAge_S#.csv
Data: flu table: SDI
*/

/* all ages ILI */
SELECT season_cdcwk.SMALL_SEAS_NUM, flu.WEEK, season_cdcwk.WK_NUM, flu.PATIENT_ZIP3, flu.AGEGROUP, flu.ILI_m from flu 
LEFT JOIN season_cdcwk USING (WEEK)
WHERE flu.SERVICE_PLACE = "TOTAL" and flu.PATIENT_ZIP3 <> "TOT" and flu.POPSTAT > 0 and season_cdcwk.SMALL_SEAS_NUM = 7
GROUP BY season_cdcwk.SMALL_SEAS_NUM, flu.WEEK, flu.PATIENT_ZIP3, flu.AGEGROUP
;


