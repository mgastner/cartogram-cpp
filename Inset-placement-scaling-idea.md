# Placement and Scaling of Insets in Output Image

## Step 1: Rescaling the Insets according to their parameters (i.e: Population, GDP, etc) to make them proportionate to eah other

### Why:

 Cartogram Insets produced by the code will not have proportionate area to each other. For example, if Inset A and Inset B have same parameter value (i.e same population, GDP), then (in a scaled version) we expect them to have same area for obvious reason; however, Cartogram Inset A and Inset B produced by the code may not have same area. So, we need to rescale in order to make them same area.

### How:

 In order for insets to be scaled relative to each other, the ratio  of Insets area and their respective parameter values (i.e population) need to be same (if Inset A areaA/parameter value is X, Inset B areaB/parameter value needs to be also X). To achieve it, we first consider the main portion of cartogram as our base and we calculate its area/parameter_value, lets call it ratioMain. Now, we need to change the Insets' area so that their parameter_value/area equals to ratioMain. After we find Insets' area, and substract sqrt(oldArea - newArea).

### Pseudocode:

```c++
//From CSV
int mainGDP, insetA_GDP, insetB_GDP;
double mainArea, insetA_Area, insetB_Area; //can be calculated using (https://github.com/mapbox/geojson-area)

//Area values we need to obtain
double insetA_Area_Scaled, insetB_Area_Scaled;

//Calculate ratio of main's parameter and Area
double ratioMain = mainArea/mainGDP;

//Calculate new Area using the ratioMain
insetA_Area_Scaled = ratioMain * insetA_GDP;
insetA_Area_Scaled = ratioMain * insetB_GDP;

//Now we know that we need to change the values inside insetA.geojson someway so that their Area becomes equal to insetA_Area_Scaled. To change insetA.geojson's area from insetA_Area to insetA_Area_Scaled, we subtract each value inside the geosjon by sqrt(inset_Area - inset_Area_Scaled).

//Therefore, we change coordinates inside of insetA.geojson-

double x_insetA_Scaled = x_insetA - sqrt(inset_Area - inset_Area_Scaled);
double y_insetA_Scaled = y_insetA - sqrt(inset_Area - inset_Area_Scaled);

```

## Step 2: Normalize the total area of the insets to be equal to 1

### Why:

It will be useful to have the actual map and cartogram's map area to be normalized to 1. Then, we can easily place them side by side or on top of each other, and their size will look similar.

### How:

It is a simple process. We take a geojson file, compute its area (all area including all insets), and subtract sqrt(Area - 1) from all of its coordinate values. "sqrt(Area - 1)" came from our equation in Step 1: sqrt(inset_Area - inset_Area_Scaled). We wanted the scaled Area to be 1.

### Pseudocode:

```c++
//from Step 1
mainArea_Scaled, insetA_Area_Scaled, insetB_Area_Scaled; //can be calculated using (https://github.com/mapbox/geojson-area)

//Compute total area
double total = mainArea_Scaled + insetA_Area_Scaled + insetB_Area_Scaled;


//Now we change coordinates of all the insets-

// For main map
x_mainArea_Scaled = x_mainArea_Scaled - sqrt(total - 1);
y_mainArea_Scaled = y_mainArea_Scaled - sqrt(total - 1);

//For insetA
x_insetA_Scaled = x_insetA_Scaled - sqrt(total - 1);
y_insetA_Scaled = y_insetA_Scaled - sqrt(total - 1);

//For insetB
x_insetB_Scaled = x_insetB_Scaled - sqrt(total - 1);
y_insetB_Scaled = y_insetB_Scaled - sqrt(total - 1);
```

## Step 3: Shift each inset so that the middle of its bounding box is at coordinate (0, 0).

### Why:

After doing so, it will be easier to position different cartogram's on top of each other. Also, this step is crucial if we want to position our insets wherever we want.

### How:

For x values, we take the average of bounding box values, xMin, xMax and subtract the average value from each x coordinate.

For y values we do the same, this time we take the average of bounding box values, yMin, yMax and subtract the average value from each y coordinate.

### Pseudocode:
```c++
//Lets consider only one inset A, even though we need to do this for all the insets and main.

//Inset A's bounding box values;
double xMin, xMax, yMin, yMax;

//We take the average of x and y values-
double xAverage = (xMin + xMax)/2;
double yAverage = (yMin + yMax)/2;

//Now, we change the coordinates of the geojson this way-
double x_New = x - xAverage;
double y_New = y - yAverage;
```

## Step 4: Find the dimensions of the combined bounding box for each inset of different parameters (i.e Population, GDP) from the normalised coordinates described above.

### Why:

Inset A bbox of Population parameter might be bigger than the Inset A bbox of GDP parameter. This won't be ideal for us if we try to reposition the insets later. Without this step, there is a possibility of overlaps.

### How:

For Inset A bboxs of population and GDP, we take the smallest bbox values that would surround two bbox's completely. We take the lesser xMin, yMin value our new xMin, yMin and bigger xMax, yMax values as our new xMax, yMax


### Pseudocode:
```c++
//Lets consider only one inset A, we have it's bounding box values for population, GDP parameter.

//Inset A's Population bounding box values;
double xMin_Pop, xMax_Pop, yMin_Pop, yMax_Pop;

//Inset A's GDP bounding box values;
double xMin_GDP, xMax_GDP, yMin_GDP, yMax_GDP;

//Calculate new bbox values for each

double xMin_both = min(xMin_Pop, xMin,GDP);
double yMin_both = min(yMin_Pop, yMin,GDP);
double xMax_both = min(xMax_Pop, xMax,GDP);
double yMax_both = min(yMax_Pop, yMax,GDP);
```

## Step 5: Add some padding to get the dimensions of the frame for each inset.

### Why:

Without it, after placing the insets, the insets might touch each other, creating the false illusion that they are geographically together.

### How:

If we want to add X amount of padding, we reduce our xMin, yMin by X, and we increase xMax, yMax by X.


### Pseudocode:
```c++
//Lets consider only one inset A, we have it's combined bounding box values from Step 4

//Inset A's combined bounding box values;
double xMin_both, xMax_both, yMin_both, yMax_both;

//Imagine we want to add P Padding;

//Calculate new bbox values with padding for each

xMin_both = xMin_both - P;
yMin_both = yMin_both - P;
xMax_both = xMax_both + P;
yMax_both = yMax_both + P;
```

## Step 6: Shift the frames to the positions specified by the users ("C", "T", "R" etc.)

### Why:

This step is necessary to shift the insets to user specified places.

### How:

We need to change the x, y values, and bbox values of the geojson considering how the user want us to place the insets (details in the pseudocode).

### Pseudocode:
```c++
//Lets consider only one inset A.

//Inset A's combined bounding box values;
double xMin_A, xMax_A, yMin_A, yMax_A;

//Inset A geojson's x, y values;
double x_InsetA, y_InsetB;

//We have the main Inset's (the one that will stay in the center) bounding box values

double xMin_main, xMax_main, yMin_main, yMax_main;

//To shift Inset A, We change the x, y values by-

//To shift the inset to the "T" (top)

x_InsetA = xInsetA;
y_InsetB = y_InsetB + (yMax_main + yMax_A)/2;

//To shift the inset to the "R" (right)

x_InsetA = xInsetA + (xMax_main + xMax_A)/2;
y_InsetB = y_InsetB;

//To shift the inset to the "B" (bottom)

x_InsetA = xInsetA;
y_InsetB = y_InsetB - (yMax_main + yMax_A)/2;

//To shift the inset to the "L" (left)

x_InsetA = xInsetA - (xMax_main + xMax_A)/2;
y_InsetB = y_InsetB;
```








