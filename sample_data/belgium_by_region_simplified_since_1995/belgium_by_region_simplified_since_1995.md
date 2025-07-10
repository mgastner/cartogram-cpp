# Sources

## belgium_by_region_simplified_since_1995.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020)
geoBoundaries: A global database of political administrative boundaries.
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866.
Downloaded from: https://github.com/wmgeolab/geoBoundaries/blob/main/releaseData/gbOpen/BEL/ADM1/geoBoundaries-BEL-ADM1_simplified.geojson on 4 July 2022.

Simplified with the following command using `mapshaper`'s command-line interface:
```bash
mapshaper ../belgium_by_region_since_1995/belgium_by_region_since_1995.geojson -simplify dp 2.5% -o belgium_by_region_simplified_since_1995.geojson
```

## simplified_belgium_population_2022.csv
National Register of Belgium. Downloaded 30 June 2022 from https://statbel.fgov.be/en/themes/population/structure-population.

Same as `../belgium_by_region_since_1995/belgium_population_2022.csv`, but with renamed to ensure different output filename.





