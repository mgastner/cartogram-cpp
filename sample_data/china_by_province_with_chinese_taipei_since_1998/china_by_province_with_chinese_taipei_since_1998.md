# Sources

## china_by_province_with_chinese_taipei_since_1998.geojson
United Nations Office for the Coordination of Humanitarian Affairs. China - Subnational Administrative Boundaries.
Downloaded from https://data.humdata.org/dataset/17a2aaa2-dea9-4a2e-8b3f-92d1bdfb850c/resource/9e67ddf9-ce26-4b7a-82b1-51e5ca0714c8/download/chn_adm_ocha_2020_shp.zip on 28 July 2022.

### Code for simplification
```
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
library(tidyverse)
temp <- tempfile()
download.file(
  "https://data.humdata.org/dataset/17a2aaa2-dea9-4a2e-8b3f-92d1bdfb850c/resource/9e67ddf9-ce26-4b7a-82b1-51e5ca0714c8/download/chn_adm_ocha_2020_shp.zip",
  temp
)
chn_unzip <- unzip(temp)
chn <- read_sf("chn_admbnda_adm1_ocha_2020.shp")
target_n_pts_in_output <- 47500
npts(chn)
chn_simp <- ms_simplify(chn, keep = target_n_pts_in_output / npts(chn))
npts(chn_simp)
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
  crs = NULL
)
unlink(chn_unzip)
  ```
  
## china_population_2020.csv
Sources: National Bureau of Statistics China. Communiqué of the Seventh National Population Census (No. 3). Downloaded from https://web.archive.org/web/20210511104847/http://www.stats.gov.cn/english/PressRelease/202105/t20210510_1817188.html on 17 July 2022.
Worldometers. Downloaded from https://www.worldometers.info/world-population/taiwan-population/, https://www.worldometers.info/world-population/china-hong-kong-sar-population/, and https://www.worldometers.info/world-population/china-macao-sar-population/ on July 17 2022.



