# Format requirements

## 10 July 2022

All data in the directory `sample_data` should be stored in subdirectories
following the naming convention:

`(geographic_region)_by_(division)_(optional_descrip)_(year[s])`

Names of geographic regions and administrative divisions should be in British
English. Please refrain from using diacritics in the names of the directories
and files. For instance, `metropolitan_france_population_by_département_2022`
should be titled `metropolitan_france_population_by_department_2022`.

Years refer to the time when the administrative divisions existed (both with
the boundaries and the names in the files stored in this directory), not the
years when the data were obtained.
It is encouraged to state the `year[s]` as a range if possible.
Examples:

- `world_by_country_from_natural_earth_since_2018`

  Explanation:
  Country borders are subject to many territorial disputes.
  Therefore, the description `from_natural_earth` is needed to indicate the
  data source.
  The go-cart.io project does not adopt a political stance on the borders
  shown.
  The latest change in the data was in 2018, when Eswatini adopted its current
  name.
  Because the data are up-to-date at the time of writing, the year should be
  stated as `since_2018`.

- `india_by_state_2000_to_2014`

  Explanation:
  Data in this directory refer to Chhattisgarh, Uttaranchal and Jharkhand (all
  formed in 2000) as separate states.
  However, Telangana (formed in 2014) does not appear as a separate state.

- `germany_by_electoral_district_2017`

  Explanation:
  Members of the German federal parliament (Bundestag) are elected every
  four years.
  Electoral districts changed in 2017 compared to 2013, and they changed again
  in 2021.
  Therefore, the directory name contains a single year (`2017`) instead of a
  range of years.

In each subdirectory, there should be:

- exactly one **GeoJSON file**, which should usually contain maximally 50,000
  points.
  The number of points can be increased if the specific example requires it
  (e.g. because there are more than 200 geographic divisions, such as counties
  or census blocks).
  The coordinates should be longitude/latitude.
  The name of the GeoJSON file should be:

  `(geographic_region)_by_(division)_(optional_descrip)_(year[s]).geojson`
  This file should thus have the same name as its parent directory.

- at least one **CSV file** for the population.
  CSV files for other statistics (e.g. GDP or members of parliament) can be
  optionally added.

  The minimal format of the CSV file should be:

  | shapeName        | Population (people) |
  | :--------------- | :------------------ |
  | Brussels-Capital | 1208542             |
  | Flanders         | 6589069             |
  | Wallonia         | 3633795             |

  - **First column** — Contains the names of the divisions (e.g., _Wallonia_).
    The values must exactly match those in the corresponding GeoJSON file.
    The column header name can be anything; however, for compatibility with other repositories (such as `cartogram-web`), it should be set to match a property key in the GeoJSON. For this example, the GeoJSON should contain a `shapeName` property with the values _Brussels-Capital_, _Flanders_, and _Wallonia_.

  - **Second column** — Contains the statistical data used as the target area for the cartogram.
    The column header should include the name and unit of the statistic (e.g., `Population (people)`).
    The unit must be enclosed in parentheses, but it may be left blank (e.g., `Population`).

  Optional (case-sensitive) columns in the CSV file:

  - **`Label` column** — Contains abbreviations for the divisions. While not necessary for cartogram generation, we **strongly recommend** including it for visualization at https://go-cart.io.

  - **`DisplayName` column** — Contains a friendly name for the divisions. This is useful if the first column uses non-friendly identifiers, such as ISO-3 country codes, but you prefer to display readable country names in the visualization at https://go-cart.io.

  There are additional examples available in the directory `sample_data`.

  CSV file names should follow the pattern:

  `(geographic_region)_(statistic)_(optional_unit)_by_(division)_(year[s]).csv`

  Years refer to the time applicable to the data, not the years when the data
  were obtained.
  As long as the GeoJSON and the CSV files refer to the same divisions, the
  `year[s]` in the file names may differ.

  Examples:

  - `world_population_by_country_2018.csv`
  - `india_agricultural_production_in_inr_by_state_2012.csv`
  - `germany_votes_for_green_party_by_electoral_district_2017.csv`

- exactly one **Markdown file** with the name:

  `(geographic_region)_by_(division)_(optional_descrip)_(year[s]).md`

  The file should include:

  - a full bibliographic reference to the data source of the geographic
    boundaries.
  - the lines of code used to convert the geographic boundaries from the
    data source to the GeoJSON in the corresponding directory.
    The code may involve subsetting, simplification, topology repair etc.
  - instructions how to run the code (e.g. stating the programming
    language, software, version of the software, version of add-on packages
    and additional files needed).
  - a full bibliographic reference to the data in the CSV file(s).

  If data are sourced from a website, the date of the download should be
  included in the reference.
