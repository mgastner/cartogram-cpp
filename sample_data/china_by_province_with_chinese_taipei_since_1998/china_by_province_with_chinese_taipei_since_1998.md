# Sources

## china_by_province_with_chinese_taipei_since_1998.geojson
Digital Map Database of China, 2020, "Provincial Boundary", https://doi.org/10.7910/DVN/DBJ3BX, Harvard Dataverse, V1.
Downloaded from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/DBJ3BX# on 17 July 2022.

### Code for simplification
```
chn <- geojson_sf("~/china_by_province_with_chinese_taipei_since_1998/province.json")
target_n_pts_in_output <- 48500
npts(chn)
chn_simp <- ms_simplify(chn, keep = target_n_pts_in_output/npts(chn))
geojson_write(
  chn_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "china_by_province_with_chinese_taipei_since_1998.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)
  ```
  
## china_population_2021.csv
Source: National Bureau of Statistics China. Communiqué of the Seventh National Population Census (No. 3). Downloaded from http://www.stats.gov.cn/english/PressRelease/202105/t20210510_1817188.html on 17 July 2022.
Worldometers. Downloaded from https://www.worldometers.info/world-population/taiwan-population/, https://www.worldometers.info/world-population/china-hong-kong-sar-population/, and https://www.worldometers.info/world-population/china-macao-sar-population/ on July 17 2022.



