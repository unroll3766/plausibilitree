# PlausibilTree

PlausibiliTree is a modular assortment of software that generates realistic and labelled point clouds of forested scenes.

## Review notes

1. utm needs to be added to requirements.txt  
2. the normalization step happens in the DL side, it’s out of scope for this project  
3. Is basis-stem still needed in quickstart instructions  
4. Please note the updated BRANCH\_ORDER spec, \-1 is used to refer to non branch points  
5. Branch order behavior is not correct. When a branch splits, one of the splits should carry the parent order. If a 0th order branch splits, the larger one needs to stay a 0th order branch while the smaller one becomes a 1st order branch.  
6. At least 2 shrubs were generated with 0 points. They have a record in individual shrubs. Diagnose why this happened.  
7. In individual\_shrubs/textual/shrub\_info.jsons:  
   1. dbh\_m can be removed entirely  
   2. null values in bb, voxel, and neighborhood can occur (likely related to 6\)  
8. in individual\_shrubs/deep\_learning/volume/labels/shrub\_volume.jsons:  
   1. Again, deep learning is out of scope. This file is not needed at all  
   2. 0 values in n\_points likely related to 6\.  
9. in individual\_trees/textual/generation.jsons:  
   1. realised values should include DBH in addition to base\_radius\_m  
10. in individual\_trees/textual/info.jsons:  
    1. convex\_hull volume is not needed (it’s populated as null anyway)  
11. Note the inclusion of a new helper module to view tree skeletons (full implementation TBD)

## Comments

For consistency, all geospatial data is ingested and exported as EPSG:4326 (lat/lon) and EPSG:5773 (elevation). Local coordinates are always measured in meters, as are parameters that operate in local space (buffers, cell sizes, etc.). Internally, Plausibilitree uses the utm python library to convert all geospatial data into a meters based projection during processing. [https\://github.com/Turbo87/utm](https://github.com/Turbo87/utm)   
We modified SimpleSynthTree to use as the tree generation module: [https\://github.com/mitchbryson/SimpleSynthTree](https://github.com/mitchbryson/SimpleSynthTree) 

## Inputs

The main input is a tree list table that represents a given plot.   
Required rows: local position (x,y coordinate pair in meters), diameter at breast height (trunk diameter at 1.37 m above ground), and tree height (highest point measured on the tree in meters), tree code (e.g. species, functional group)  
Optional rows: lat/lon (in EPSG:4326), elevation (EPSG:5773)

A tree code \-\> simulator parameter dictionary. These are software specific and contain the desired parameters for a tree with a given DBH, HT, and tree code.

User selected generator module (currently just SST)

Topography dict (DTM path optional)

Vegetation parameter dict

Point cloud settings (max density)

## Process

* Starting with the first tree: generate a tree using the appropriate parameters given the DBH, HT, and tree code from the tree parameter dict. Place it at the provided x,y with a random rotation \[0,360) in steps of ROT\_STEPS and a random reflection. Continue until all trees are placed  
* Use the tree skeletons to identify intersections with other trees. For each tree, in ascending order of UID:  
  * When parts from two different trees intersect, points from one of the trees are flagged as trimmed=TRUE. Rules for trim selection  
    * When points are flagged for trimming, find the nearest joint before the collision. Any points connected to that joint are flagged as trimmed=TRUE.  
    * If two branches intersect, the one with the shortest branch length (measured as maximum length of connected segments from an intersecting point back to trunk) is trimmed. In the case of a tie, randomly select a branch to cut  
    * If a branch intersects with a trunk, the branch is trimmed.  
    * If two trunks intersect, the trunk with the shortest skeletal length is trimmed. In the case of a tie, select the trunk associated with the higher Tree ID to cut.  
  * Any points flagged from trimming are removed from the point cloud.  
* Next, the program generates topography. If a DTM was provided, trim it to plot boundary and random sample from it to generate ground points. If a DTM was not provided, use the random terrain parameters to generate a DTM and sample from it. If none was provided, random sample ground points from a uniform plane.  
  * If multiple options are provided (filepath, random parameters, and/or none), resolve it in this way: None if None was included, else choose filepath if filepath was provided. Return a warning describing this choice.   
  * If none was provided and sampling density was not provided, return a warning and sample at a default density.  
  * If the topography dict is otherwise malformed or not provided, return a warning and operate as if None was provided.  
* TBD implementation: adding boulders and rocks to the plot in unoccupied areas.  
  * User specifies the number and size range of rocks and boulders. There also needs to be some kind of parameter to account for how irregular a boulder might be.  
  * Boulders and rocks are placed the same as trees. The lowest part of the boulder is placed at the lowest part of the ground, thus any boulders and rocks “sink” into the ground on the upslope.  
* Next, the program generates vegetation using the vegetation parameter dict.   
  * Shrubs are populated first. Each shrub is randomly generated (sampled from shrub parameters). The shrub is placed randomly within the constraints of the shrub pattern and tree proximity. Any point belonging to a shrub that intersects with a tree point is removed.  
  * Grass is populated next. Any point belonging to grass that intersects with a tree point or a shrub point is removed.  
* Finally, the program populates downed woody debris. TBD how to implement this but::  
  * Downed logs of varying radius, length  
  * Tree stumps of varying radius and height  
  * Stumps need to have a minimum distance from other trees  
  * Downed logs need to have a minimum and maximum distance from other trees  
  * Number of stumps, number of logs

## Outputs

Plot level point cloud (.laz) file with the following extra dimensions:

* Semantic code (trunk, branch, foliage, shrub, grass, ground)  
* Tree UID  
* Tree code   
* Shrub UID  
* Height above ground (height of every vegetation point above the DTM)

Plot level  
JSON: number of trees, number of shrubs, volume of trunks, branches, foliage, shrubs, and grass, bounding box  
DTM (if randomly generated)  
Tree level   
JSON: tree ID, position, volume of trunk, branch, foliage, shrub, grass, bounding box, tree DBH, tree height

# Helper tools 

## FIA inventory \-\> tree list table

The FIA helper leverages widely available data (FIA) to rapidly create plausible tree lists. Trees do not come with individual coordinates assigned, and variable plot size means trees are not recorded exhaustively over the plot. This helper works at the scale of one unit of analysis (i.e. a plot, or a plot+subplot combination)

### Inputs

various FIA tables for a state  
user selects plot shape (rectangular, circular) and dimensions (length/width, radius)  
user selects expansion factor modifier (modifier is multiplied by FIA recorded TPA). 0 means that only documented trees are generated, values between 0 and 1 (exclusive) are a global reduction to FIA TPA, greater than 1 means more trees will be generated than the given TPA.  
user selects minimum distance between trees (negative numbers means trees can overlap, 0 means trees can’t touch, and positive numbers mean trees must have some buffer from each other)   
user selects maximum number of tries (FIA helper will attempt to place a tree up to maximum number of tries)  
user selects edge effects, which ranges from 0.0 to 1.0. The proportion of a tree trunk’s area contained outside the plot must be less than or equal to edge effects. For instance, at edge effect 1.0, up to 100% of the tree trunk may fall outside the plot, while at edge effect 0.0, no part of the tree trunk may fall outside the plot. Note: while FIA helper will never place a tree 100% outside the plot, giving users the full 0 to 1 range allows them to ignore the geometry of edge effects.

### Process

* User can select to generate lists using plots or subplots  
* FIA helper extracts the geographic coordinate from plot center  
* Remove any trees with DECAYCD \>= 1 (snags out of scope)  
* Remove any trees with missing data in the following columns  
  * DIA, HT, TPA\_UNADJ, CCLCD, selected functional group column (SPCD, SPGRPCD, etc.)  
* Convert the DIA column to meters  
* Convert the HT column to meters  
* Convert the TPA\_UNADJT column to square meters, multiply it by the user provided modifier, and round it up to the nearest integer.  
* FIA helper uses the expansion factor (and the user provided modifier) to populate the number of trees expected for the user supplied plot dimensions.  
* For each tree, check TPA\_UNADJT  
  * If N \> 1, duplicate that row N times  
  * If N \= 1, do nothing  
  * If N \= 0, remove the tree from the dataframe  
* FIA helper randomly places the trees within the user supplied plot dimensions.  
  * Sort the trees to be placed in descending order of height, with DBH as a tie breaker. If two trees have the same height and DBH, randomly select one to be first.  
  * Randomly select an x,y position contained within the user supplied plot dimensions. That will be the coordinates (centroid) of the first tree.  
  * Calculate the trunk area of the placed tree by creating a circle with diameter (DBH+minimum distance between trees) placed on the centroid  
  * Check if the trunk area overlaps with any previously placed trees. If so, do not place the tree and increment the current number of tries. If the current number of tries is less than max number of tries, select a new random position and try again. Otherwise, do not place the tree.  
  * Check if the proportion of the trunk area outside the plot is less than or equal to edge effects. If so, do not place the tree and increment the current number of tries. If the current number of tries is less than max number of tries, select a new random position and try again. Otherwise, do not place the tree.  
  * If the tree is successfully placed, remove it from the list, generate a unique ascending order tree UID and append its position and FIA information to the output tree list.  
  * Repeat until the tree placement list is empty.

### Outputs:

a .csv of all successfully placed trees with the following:  
Required rows: local position (x,y coordinate pair in meters), diameter at breast height (trunk diameter at 1.37 m above ground), and tree height (highest point measured on the tree in meters), tree code (e.g. species, functional group)  
Optional rows: lat/lon (local position transformed into EPSG:4326 using local position and geographic coordinates of plot center), elevation (meters above sea level)  
a JSON that contains the parameters used in FIA helper, the number of trees successfully placed, the number of trees that could not be placed, and the proportion.

## Helper: Tree Parameters (SimpleSynthTree)

This is a software tool to help users establish parameters that produce realistic trees for a given tree code.

### Inputs:

A tree list .CSV with the following required rows: DBH, HT, tree code  
Optional input: a JSON parameter dict.

### Process

User provides a tree list and/or a previous parameter dict  
User can select a tree from the tree list, or create one by hand (e.g. manual input of DBH, HT, tree code)  
User can edit the parameter dict values if provided. If one wasn’t provided, program generates a table with a row for each unique tree code in the tree list.   
User can click a generate button that creates the selected tree in SST using the current parameter dict. The program has a 3D viewer pane where the user can inspect the tree.  
User can save and export the modified parameter dict as a JSON.

### Outputs

A JSON parameter dict.

## Helper: Vegetation Parameters (all modules)

### Inputs:

Trimmed tree plot  
Vegetation parameter dict (optional)

### Process:

Program has a 3D viewer pane for the original trimmed tree plot  
User can edit the vegetation dict values if provided. If one wasn’t provided, program loads a dict with default values.  
User can select a generate button that adds in the vegetation points to the chosen cloud. This is displayed in the 3D viewer pane. User should be able to toggle between the original and vegetation added plot.  
User can save and export the current vegetation dict.

### Output:

A JSON parameter dict.

## Helper: HELIOS++ export

Uses HELIOS++ to simulate scanning effects

### Input

Point cloud (voxelized) or polygons   
Scanner parameters

### Process

### Output

A .laz file generated by HELIOS++

## Helper: Tree skeleton viewer

A 3D GUI to load a tree skeleton JSON

### Input

skeleton.json (raw or trimmed)

### Process

User should be able to toggle various attributes:

* Radius \- displays the width of each branch when activated. Shows polylines when deactivated.  
* Branch ID \- toggle a display on or off.  
* Branch angle \- toggle a display on or off.  
* Branch order \- colorize skeleton/cylinders by branch order.

User should be able to zoom in/out, rotate the skeleton, etc.  
User should be able to render the current view to a .png or .jpeg

### Output

TBD

## SimpleSynthTree modifications

### Tree generation modifications

* Added a translation factor so tree tops can grow towards any point  
* Parameter to generate multiple branches   
* Sampling density change. User provides the number of points used for a height slice (256, 512, etc.) and the program uses the height of the tree and its DBH to calculate the total number of points (more DBH \= more points)

### Skeleton export modification

* Added an option to export a JSON that describes how branches are connected to one another and their volume

### .LAS update modification

* Updated to 1.4 format  
* Added additional point columns: tree\_ID, semantic, branch\_ID, branch\_order  
* Add support for .laz export

# Variables

PLOT\_SHAPE \= string. Dictates whether the plot is circular or rectangular. Must be one of “CIRCLE” or “RECTANGLE”, defaults to “CIRCLE”  
(PLOT\_LENGTH, PLOT\_WIDTH) \= tuple of two floats. length and width (in meters) of a rectangular plot, defaults to (30.0, 30.0). Only valid when PLOT\_SHAPE \== “RECTANGLE”  
PLOT\_RADIUS \= float. radius of a circular plot (in meters), defaults to 7.315. Only valid when PLOT\_SHAPE \== “CIRCLE”  
MIN\_DBH \= float. Minimum size for a tree to be placed (in meters). Defaults to 0.127.  
DBH\_HT \= float. Height at which diameter is measured for DBH (in meters). Defaults to 1.37  
FOLIAGE\_VOXEL \= float. Size used for voxel downsample to calculate foliage volume (in meters). Defaults to 0.1  
GRASS\_VOXEL \= float. Size used for voxel downsample to calculate grass volume (in meters). Defaults to 0.1  
DTM\_RESOLUTION \= float. Resolution used to generate the DTM. Defaults to 0.25  
CHM\_RESOLUTION \= float. Resolution used to generate the canopy height map. Defaults to 0.05  
STEM\_RESOLUTION \= float. Resolution used to generate the stem map. Defaults to 0.05  
STEM\_SECTION \= float. The buffer applied to DBH\_HT when generating a stem map. Defaults to 0.25.  
ROT\_STEPS \= int. Number of steps (in degrees) to take when rotating a tree, defaults to 30\.

# Data Structures

## Tree List (Input to Plausibilitree)

Required attributes:

| Attribute | Data type | Data description | Range | Notes |
| :---- | :---- | :---- | :---- | :---- |
| TREE\_UID | string | A unique ID for each tree in the list |  | ex: 0001, 0002 |
| X | float | local x coordinate of centroid of tree, in meters, measured from lower left of plot | \[0, PLOT\_LENGTH\] (Rectangular plot) or \[0, 2\*PLOT\_RADIUS\] (circular plot) |  |
| Y | float | local y coordinate of centroid of tree, in meters, measured from lower left of plot | \[0, PLOT\_WIDTH\] (Rectangular plot) or \[0, 2\*PLOT\_RADIUS\] (circular plot) |  |
| Z | float | local z coordinate of centroid of tree, in meters, measured from base of tree | \[0, inf) | typically will be 0 for all trees |
| DBH | float | diameter of the tree trunk at breast height (DBH\_HT), in meters | (MIN\_SIZE, PLOT\_WIDTH\] (Rectangular plot) or (0, 2\*PLOT\_RADIUS\] (circular plot) |  |
| HT | float | height of the tree above ground, measured in meters | (0, 120\) | Tallest tree in the world is \~116 m tall |
| TREE\_CODE | string | a valid code that describes the tree’s vegetation type (species code, functional group code, etc.) |  | See FIA codes for examples |

Optional attributes

| Attribute | Data type | Data description | Range | Notes |
| :---- | :---- | :---- | :---- | :---- |
| LON | float | longitude of the tree’s centroid, measured in decimal degrees, in EPSG:4326 | \[-180.0,180.0\] |  |
| LAT | float | latitude of the tree’s centroid, measured in decimal degrees, in EPSG:4326 | \[-90.0,90.0\] |  |
| ELEV | float | elevation of the base of the tree’s centroid, measured in meters, in EPSG:5773  | \[-500, 5000\] | Range estimated with AI |

## Plot level outputs of Plausibilitree

### Polygon (.ply)

Plot level mesh with the following elements:

* Semantic code (1=trunk, 2=branch, 3=foliage, 8=shrub, 9=grass, 0=no semantic class (ground, noise, etc.)  
* Tree UID (string)  
* Tree code (string)  
* Shrub UID (string)

### Point Cloud (.laz)

Plot level point cloud with the following extra dimensions:

* Semantic code (1=trunk, 2=branch, 3=foliage, 8=shrub, 9=grass, 0=ground)  
* Tree UID (string)  
* Tree code (string)  
* Shrub UID (string)  
* Height above ground (height of every point above the DTM) (float)

### Plot info (.JSON)

Required attributes:

| Attribute | Data type | Data description | Lower range | Upper range | Notes |
| :---- | :---- | :---- | :---- | :---- | :---- |
| NUM\_TREES | int | Number of trees contained in the plot | \>= 0 | inf |  |
| BASAL\_AREA | float | Basal area of trees contained in plot | \>= 0 | inf | measured in square meters (use 0 if no trees) |
| NUM\_SHRUBS | int | Number of shrubs contained in the plot | \>= 0 | inf |  |
| VOL\_TREES | float | Sum of the volume of every tree | \>= 0 | PLOT\_WIDTH\*PLOT\_LENGTH\*MAX(TREE\_LIST\$HT) | Volume comes from individual tree .JSON files |
| VOL\_SHRUB | float | Sum of the volume of every tree | \>= 0 | PLOT\_WIDTH\*PLOT\_LENGTH\* MAX(SHRUB\_LIST\$HT) | Volume comes from individual shrub .JSON files.  |
| VOL\_GRASS | float | Sum of the volume of all points labelled grass | \>= 0 | PLOT\_WIDTH\*PLOT\_LENGTH\*MAX(TREE\_LIST\$HT) | Volume calculated from grass points voxelized at GRASS\_VOXEL |
| VOL\_TRUNK | float | Sum of the volume of all portions labelled trunk | \>= 0 | PLOT\_WIDTH\*PLOT\_LENGTH\*MAX(TREE\_LIST\$HT) | Volume calculated from individual tree .JSON files |
| VOL\_BRANCH | float | Sum of the volume of all portions labelled branch | \>= 0 | PLOT\_WIDTH\*PLOT\_LENGTH\*MAX(TREE\_LIST\$HT) | Volume calculated from individual tree .JSON files |
| VOL\_FOLIAGE | float | Sum of the volume of all points labelled foliage | \>= 0 | PLOT\_WIDTH\*PLOT\_LENGTH\*MAX(TREE\_LIST\$HT) | Volume calculated from foliage points voxelized at FOLIAGE\_VOXEL |
| PLOT\_SHAPE | string | one of circle or rectangle |  |  |  |
| PLOT\_AREA | float | area | \> 0 | inf |  |
| BOUNDING\_BOX | 6 floats | 3D bounding box of point cloud |  |  | If global coordinates are available, use EPSG:4326 for x,y and EPSG:5773 for z. Otherwise, use local coordinates in meters. |
| CHM\_PATH | string | filepath to the generated CHM |  |  |  |
| STEM\_PATH | string | filepath to the generated stem map |  |  |  |
| PLAUSIBILITREE\_CONFIGURAITION | string | filepath to .JSON recording settings used |  |  |  |

Optional attributes:

| Attribute | Data type | Data description | Lower range | Upper range | Notes |
| :---- | :---- | :---- | :---- | :---- | :---- |
| PLOT\_CENTER | (float, float, float) | Lon, lat, and elevation of plot center according to EPSG:4326 (lon/lat), and EPSG:5773 (elev) | \-180, \-90, \-500 | 180, 90, 5000 | Range orders in (lon, lat, elev). Elevation range estimated with AI |
| DTM\_PATH | string | filepath to generated DTM  |  |  |  |

### Canopy height map (.geoTIFF)

A canopy height map in EPSG:4326 with a resolution of CHM\_RESOLUTION and 2 bands. The first band contains the max height above ground for all vegetation points within the cell (0 if no vegetation). The second band is a mask that contains 1 if foliage points are within the cell, or 0 if no foliage points are present.

### Stem map (.geoTIFF)

A stem map in EPSG:4326 with a resolution of STEM\_RESOLUTION and 2 bands. The first band contains the count of all points with a height above ground greater than or equal to DBH\_HT-STEM\_BUFFER and less than or equal to DBH\_HT+STEM\_BUFFER. The second band contains 1 for all tree points with a height above ground greater than or equal to DBH\_HT-STEM\_BUFFER and less than or equal to DBH\_HT+STEM\_BUFFER and a 0 anywhere else.

### DTM (.geoTIFF)

A digital terrain map in EPSG:4326 with a resolution of DTM\_RESOLUTION. Only required if the user generated a DTM. 

## Tree level outputs of Plausibilitree

### Tree info (.JSON)

Required attributes:

| Attribute | Data type | Data description | Lower range | Upper range | Notes |
| :---- | :---- | :---- | :---- | :---- | :---- |
| TREE\_UID | string | A unique ID for each tree in the list |  |  | Should correspond to the tree UID from tree list |
| TREE\_DBH | float | diameter of the tree trunk at breast height (DBH\_HT), in meters | \<= 0 | PLOT\_WIDTH (Rectangular plot) or 2\*PLOT\_RADIUS (circular plot) | May differ from DBH in tree list |
| TREE\_HT | float | height of the tree above ground, measured in meters | \<= 0 | 120 | May differ from HT in tree list. Tallest tree in the world is \~116 m tall |
| VOL\_TREE | float | Sum of the volume of trunk, branch, and foliage | \<= 0 |  |  |
| VOL\_TRUNK | float | Sum of the volume of all portions labelled trunk | \<= 0 |  | Volume of cylinders corresponding to trunk |
| VOL\_BRANCH | float | Sum of the volume of all portions labelled branch | \<= 0 |  | Volume of cylinders corresponding to branch |
| VOL\_FOLIAGE | float | Sum of the volume of all points labelled foliage | \<= 0 |  | Volume calculated from that tree’s foliage points voxelized at FOLIAGE\_VOXEL |
| BOUNDING\_BOX | 6 floats | 3D bounding box of point cloud |  |  | If global coordinates are available, use EPSG:4326 for x,y and EPSG:5773 for z. Otherwise, use local coordinates in meters. |
| TREE\_X | float | local x coordinate of centroid of tree, in meters, measured from lower left of plot | 0 | PLOT\_LENGTH (Rectangular plot) or 2\*PLOT\_RADIUS (circular plot) | Should be same as tree list |
| TREE\_Y | float | local y coordinate of centroid of tree, in meters, measured from lower left of plot | 0 | PLOT\_LENGTH (Rectangular plot) or 2\*PLOT\_RADIUS (circular plot) | Should be same as tree list |
| TREE\_Z | float | local z coordinate of centroid of tree, in meters, measured from base of tree | \[0, inf) | \[0, inf) | This can, and often will, vary from the tree list z |

Optional attributes:

| Attribute | Data type | Data description | Lower range | Upper range | Notes |
| :---- | :---- | :---- | :---- | :---- | :---- |
| TREE\_CENTER\_REAL | (float, float, float) | Lon, lat, and elevation of tree’s centroid according to EPSG:4326 (lon/lat), and EPSG:5773 (elev) | \-180, \-90, \-500 | 180, 90, 5000 | Range orders in (lon, lat, elev). Elevation range estimated with AI. Use bottom of tree for elev |

### Tree skeleton (.JSON)

A .JSON where each element represents a portion of the generated tree.

| Attribute | Data type | Data description | Lower range | Upper range | Notes |
| :---- | :---- | :---- | :---- | :---- | :---- |
| BRANCH\_UID | string | A unique ID for each element in the tree |  |  |  |
| BRANCH\_ORDER | int | A 0th order branch is the trunk, a 1st order branch is connected to the trunk, a 2nd order branch is connected to a 1st order banch, and so on. \-1 indicates the point does not belong to a branch | \-1 | 5 | 5 is a soft upper limit. |
| BRANCH\_LENGTH | float | The length of the element | \> 0 |  |  |
| BRANCH\_RADIUS | float | The base radius of the element | \> 0 |  |  |
| BRANCH\_VOLUME | float | The volume of a cylinder fitted calculated with the element’s length and radius | \> 0 |  |  |
| BRANCH\_ANGLE | (float, float, float) | The angle of the branch w.r.t. x, y, z axis in degrees | \>= 0 | \< 360 |  |
| SOURCE | string | UID of the element this element connects to |  |  | Every element should only have one source |
| DESTINATION | string | UID(s) of the elements connected to this element |  |  | In the case of a split branch, there may be multiple destinations |

## Tree Parameter Dict

These parameters are used to convert DBH, HT for a given tree code into SimpleSynthTree parameters

* a0 and a1 (crown width)  
* taper function (base diameter)  
* DBH min/max (these are bounds on the tree’s DBH)  
* HT min/max (these are bounds on the tree’s HT)

These parameters are required for SimpleSynthTree

* HT\_LO, HT\_HI  
* BASE\_DIA\_LO, BASE\_DIA\_HI  
* SPLIT\_HT\_LO  
* SPLIT\_HT\_HI  
* SPLIT\_PROB  
* NUM\_BRANCH\_LO  
* NUM\_BRANCH\_HI  
* MIN\_CAN\_HT\_LO  
* MIN\_CAN\_HT\_HI  
* MAX\_CAN\_W  
* MAX\_CAN\_W\_HT\_LO  
* MAX\_CAN\_W\_HT\_HI  
* TOP\_DIST  
* MID\_DIST  
* FP\_RAND\_HT  
* SPR\_LO  
* SPR\_HI

## Vegetation Dict

* Shrubs  
  * Shrub cover  
  * Min/max horizontal size  
  * Min/max shrub height  
  * Shape complexity  
  * Density/spatial pattern  
  * Tree proximity  
* Grass  
  * Height  
  * Density  
  * Spatial pattern  
  * Tree proximity 

## Topography Dict

* One of  
  * DTM filepath   
  * Parameters for random terrain generation  
  * None  
* Ground sampling density (DTM\_RESOLUTION)

## 
