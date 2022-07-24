# Sources

## china_by_province_with_chinese_taipei_since_1998.geojson
Global Administrative Areas (2012). GADM database of Global Administrative Areas, version 2.0. [online] URL: www.gadm.org.
Downloaded from https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_CHN_shp.zip, https://earthworks.stanford.edu/catalog/stanford-yq344jw1449, https://geodata.lib.utexas.edu/catalog/stanford-hh544xg3454, and https://github.com/wmgeolab/geoBoundaries/blob/main/releaseData/gbOpen/TWN/ADM1/geoBoundaries-TWN-ADM1_simplified.geojson on 24 July 2022.

### Code for simplification
```
chn <- geojson_sf("~/china_by_province_with_chinese_taipei_since_1998/gadm36_CHN_1.json")
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
  
## china_population_2020.csv
Sources: National Bureau of Statistics China. Communiqué of the Seventh National Population Census (No. 3). Downloaded from https://web.archive.org/web/20210511104847/http://www.stats.gov.cn/english/PressRelease/202105/t20210510_1817188.html on 17 July 2022.
Worldometers. Downloaded from https://www.worldometers.info/world-population/taiwan-population/, https://www.worldometers.info/world-population/china-hong-kong-sar-population/, and https://www.worldometers.info/world-population/china-macao-sar-population/ on July 17 2022.



