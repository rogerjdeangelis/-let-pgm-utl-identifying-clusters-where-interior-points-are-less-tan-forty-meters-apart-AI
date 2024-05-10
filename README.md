# -let-pgm-utl-identifying-clusters-where-interior-points-are-less-tan-forty-meters-apart-AI
Identifying clusters where interior points are less tan 40 meters apart AI 
    %let pgm=utl-identifying-clusters-where-interior-points-are-less-tan-forty-meters-apart-AI;

    Identifying clusters where interior points are less tan 40 meters apart AI


    Output cluster plot
    https://tinyurl.com/3ys2efhj
    https://github.com/rogerjdeangelis/utl-identifying-clusters-where-interior-points-are-less-tan-forty-meters-apart-AI/blob/main/clusters

    github
    https://tinyurl.com/2j4vhzdn
    https://github.com/rogerjdeangelis/utl-identifying-clusters-where-interior-points-are-less-tan-forty-meters-apart-AI

    inspired by
    https://goo.gl/seJSjr
    http://gis.stackexchange.com/questions/17638/how-to-cluster-spatial-data-in-r

    You need to install these

    install.packages("sp");
    install.packages("geosphere");
    install.packages("dismo");
    install.packages('rgdal', type="source")
    remotes::install_version("rgeos", version = "0.6-4")
    install.packages("https://cran.r-project.org/src/contrib/Archive/rgdal/rgdal_1.4-8.tar.gz", repos=NULL, type="source")


    Related repos

    https://github.com/rogerjdeangelis/utl-R-kmeans-clustering--of-tweets
    https://github.com/rogerjdeangelis/utl-draw-polygons-around-clusters-of-points-r-convex-hulls
    https://github.com/rogerjdeangelis/utl-kmeans-cluster-analysis-based-on-only-one-variable-in-r
    https://github.com/rogerjdeangelis/utl-optimum-number-of-clusters-elbow-plot
    https://github.com/rogerjdeangelis/utl-simple-informative-cluster-analysis-using-sas-interface-to-R
    https://github.com/rogerjdeangelis/utl_grouping_monthly_checking_account_into_two_clusters_by_year
    https://github.com/rogerjdeangelis/utl_linked_and_unlinked_clusters_of_servers_that_have_one_or_more_linkages_in_common

    /*               _     _
     _ __  _ __ ___ | |__ | | ___ _ __ ___
    | `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \
    | |_) | | | (_) | |_) | |  __/ | | | | |
    | .__/|_|  \___/|_.__/|_|\___|_| |_| |_|
    |_|
    */

    /**************************************************************************************************************************************
    /*                                              |                                   |                                                 *
    /*                     INPUT                    |           PROCESS                 |                       OUTPUT                    *
    /*                                              |                                   |                                                 *
    /*                        X                     | %utl_submit_r64('                 |                           X                     *
    /*       -1.486  -1.485      -1.483      -1.482 | library(sp);                      |           -1.485   -1.484   -1.483   -1.482     *
    /*       --+---//--+---------/-+------------+-  | library(rgdal);                   |       ------+--------+--------+--------+----+   *
    /*    Y  |                                   |  | library(geosphere);               |    Y |                                      |   *
    /*       |                                   |  | library(dismo);                   |      |  Because the diameter      ******    |   *
    /*   54.9+                                   +  | library(rgeos);                   | 54.90+  of each circle          ***    ***  +   *
    /*       |                                X  |  | library(haven);                   |      |  is 40 meters.          **        ** |   *
    /*       |                              X X  |  | have=read_sas(                    | 54.90+  The distance           *   XX     * +   *
    /*       |                                   |  |   "d:/sd1/have.sas7bdat"          |      |  beteen points          *          * |   *
    /*  54.90+                                   +  |   );                              | 54.90+  are all less tha       ***  X   *** +   *
    /*       |                                   |  | x<-have$X;                        |      |  40 meters.               ********   |   *
    /*       |                                   |  | y<-have$Y;                        | 54.90+  Uses cluster                        +   *
    /*       |                       X           |  | xy <- SpatialPointsDataFrame(     |      |  centraoid          ******           |   *
    /*     54+                                   +  |  matrix(c(x,y)                    | 54.90+                   ***    ***         +   *
    /*       |                                   |  |  ,ncol=2)                         |      |                  **        **        |   *
    /*       |                                   |  |  ,data.frame(ID=seq(1:length(x))) | 54.89+                  *    X     *        +   *
    /*       |                                   |  |  ,proj4string                     |      |                  *          *        |   *
    /*  54.89+                                   +  |  =CRS(                            | 54.89+                  ***      ***        +   *
    /*       |  X X                              |  |  "+proj=longlat                   |      |                    ********          |   *
    /*       |                                   |  |   +ellps=WGS84                    | 54.89+    ******                            +   *
    /*       |                                   |  |   +datum=WGS84"                   |      |  ***    ***                          |   *
    /*   54.8+   X                               +  |    ));                            | 54.89+ **  X XX  **                         +   *
    /*       |                                   |  |  "xy";                            |      | *-- 40m ---*                         |   *
    /*       |                  X                |  |   xy;                             | 54.89+ *    X     *                         +   *
    /*       |                                   |  | mdist <- distm(xy);               |      | ***      ***      ******             |   *
    /*  54.89+                                   +  | "mdist";                          | 54.89+   ********      ***    ***           +   *
    /*       |                                   |  | mdist;                            |      |                **        **          |   *
    /*       --+----//-+---------/-+------------+-  | hc <- hclust(as.dist(mdist)       | 54.89+                *    X     *          +   *
    /*       -1.486  -1.485  X   -1.483      -1.482 | , method="complete");             |      |                *          *          |   *
    /*                                              | hc;                               | 54.89+                ***      ***          +   *
    /*                                              | d=40;                             |      |                  ********            |   *
    /*  libname sd1 "d:/sd1";                       | xy$clust <- cutree(hc, h=d);      |      ------+--------+--------+--------+-----    *
    /*  data sd1.have;                              | "xy$clust";                       |          -1.486   -1.485   -1.484   -1.483      *
    /*   array xs[9] (-1.482156,-1.482318           | xy$clust;                         |                          X                      *
    /*   ,-1.482129, -1.482880, -1.485735           | xy@bbox[] <-                      |                                                 *
    /*   ,-1.485770, -1.485913, -1.484275           |   as.matrix(extend(extent(xy)     |                                                 *
    /*   ,-1.485866);                               |  ,0.001));                        |                                                 *
    /*   array ys[9] (54.90083, 54.90078            | xy@bbox[];                        |                                                 *
    /*  ,54.90077, 54.90011, 54.89936               | cent <- matrix(ncol=2             |                                                 *
    /*  ,54.89935, 54.89935, 54.89879               |    ,nrow=max(xy$clust));          |                                                 *
    /*  ,54.89902);                                 | for (i in 1:max(xy$clust))        |                                                 *
    /*   do i=1 to 9;                               |     cent[i,] <-                   |                                                 *
    /*     x=xs[i];                                 |          gCentroid(subset(xy      |                                                 *
    /*     y=ys[i];                                 |          ,clust == i))@coords;    |                                                 *
    /*     keep x y;                                | ci <- circles(cent                |                                                 *
    /*     output;                                  |  ,d=d, lonlat=T);                 |                                                 *
    /*   end;                                       | "cent";                           |                                                 *
    /*  run;quit;                                   | cent;                             |                                                 *
    /*                                              | ci;                               |                                                 *
    /*      X          Y                            | pdf("d:/pdf/clusters.pdf");       |                                                 *
    /*                                              | plot(ci@polygons, axes=T);        |                                                 *
    /*  -1.48216    54.9008                         | ci@polygons;                      |                                                 *
    /*  -1.48232    54.9008                         | plot(xy, col=                     |                                                 *
    /*  -1.48213    54.9008                         | rainbow(4)[factor(xy$clust)]      |                                                 *
    /*  -1.48288    54.9001                         |,add=T);                           |                                                 *
    /*  -1.48574    54.8994                         |');                                |                                                 *
    /*  -1.48577    54.8994                         |                                   |                                                 *
    /*  -1.48591    54.8994                         |                                   |                                                 *
    /*  -1.48428    54.8988                         |                                   |                                                 *
    /*  -1.48587    54.8990                         |                                   |                                                 *
    /*                                              |                                   |                                                 *
    /**************************************************************************************************************************************

    /*                   _
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */

    options validvarname=upcase;
    libname sd1 "d:/sd1";
    data sd1.have;
     array xs[9] (-1.482156, -1.482318, -1.482129, -1.482880, -1.485735, -1.485770, -1.485913, -1.484275, -1.485866);
     array ys[9] (54.90083, 54.90078, 54.90077, 54.90011, 54.89936, 54.89935, 54.89935, 54.89879, 54.89902);
     do i=1 to 9;
       x=xs[i];
       y=ys[i];
       keep x y;
       output;
     end;
    run;quit;

    options ls=64 ps=24;
    proc plot data=sd1.have(rename=y=y123456789012345678901234567890);
      plot y123456789012345678901234567890*x="X" /box;
    run;quit;
    options ls=255 ps=66;

    /***************************************************************************************************************************/
    /*                                                                                                                         */
    /*                                                        X                                                                */
    /* SD1.HAVE total obs=9                 -1.486        -1.484        -1.482                                                 */
    /*                                    Y    -+-------------+-------------+--                                                */
    /* Obs        X          Y                 |                              |                                                */
    /*                                    54.901 +                            + 54.90                                          */
    /*  1     -1.48216    54.9008              |                          XX  |                                                */
    /*  2     -1.48232    54.9008              |                              |                                                */
    /*  3     -1.48213    54.9008              |                              |                                                */
    /*  4     -1.48288    54.9001         54.900 +                      X     + 54.90                                          */
    /*  5     -1.48574    54.8994              |                              |                                                */
    /*  6     -1.48577    54.8994              |                              |                                                */
    /*  7     -1.48591    54.8994              | XX                           |                                                */
    /*  8     -1.48428    54.8988         54.899 + X                          + 54.89                                          */
    /*  9     -1.48587    54.8990              |            X                 |                                                */
    /*                                         |                              |                                                */
    /*                                         |                              |                                                */
    /*                                    54.898 +                              54.89+                                         */
    /*                                         |                              |                                                */
    /*                                         -+-------------+-------------+--                                                */
    /*                                       -1.486        -1.484        -1.482                                                */
    /*                                                        X                                                                */
    /*                                                                                                                         */
    /***************************************************************************************************************************/

    /*
     _ __  _ __ ___   ___ ___  ___ ___
    | `_ \| `__/ _ \ / __/ _ \/ __/ __|
    | |_) | | | (_) | (_|  __/\__ \__ \
    | .__/|_|  \___/ \___\___||___/___/
    |_|
    */

    %utl_submit_r64('
    library(sp);
    library(rgdal);
    library(geosphere);
    library(dismo);
    library(rgeos);
    library(haven);
    have=read_sas("d:/sd1/have.sas7bdat");
    x<-have$X;
    y<-have$Y;
    xy <- SpatialPointsDataFrame(
          matrix(c(x,y), ncol=2), data.frame(ID=seq(1:length(x))),
          proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"));
    "xy";
    xy;
    mdist <- distm(xy);
    "mdist";
    mdist;
    hc <- hclust(as.dist(mdist), method="complete");
    hc;
    d=40;
    xy$clust <- cutree(hc, h=d);
    "xy$clust";
    xy$clust;
    xy@bbox[] <- as.matrix(extend(extent(xy),0.001));
    xy@bbox[];
    cent <- matrix(ncol=2, nrow=max(xy$clust));
    for (i in 1:max(xy$clust))
        cent[i,] <- gCentroid(subset(xy, clust == i))@coords;
    ci <- circles(cent, d=d, lonlat=T);
    "cent";
    cent;
    ci;
    pdf("d:/pdf/clusters.pdf");
    plot(ci@polygons, axes=T);
    ci@polygons;
    plot(xy, col=rainbow(4)[factor(xy$clust)], add=T);
    ');

    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /* d:/pdf/clusters.pdf                                                                                                    */
    /*                                                                                                                        */
    /*                           X                                                                                            */
    /*           -1.485   -1.484   -1.483   -1.482                                                                            */
    /*       ------+--------+--------+--------+----+                                                                          */
    /*    Y |                                      |                                                                          */
    /*      |  Because the diameter      ******    |                                                                          */
    /* 54.90+  of each circle          ***    ***  +                                                                          */
    /*      |  is 40 meters.          **        ** |                                                                          */
    /* 54.90+  The distance           *   XX     * +                                                                          */
    /*      |  beteen points          *          * |                                                                          */
    /* 54.90+  are all less tha       ***  X   *** +                                                                          */
    /*      |  40 meters.               ********   |                                                                          */
    /* 54.90+  Uses cluster                        +                                                                          */
    /*      |  centraoid          ******           |                                                                          */
    /* 54.90+                   ***    ***         +                                                                          */
    /*      |                  **        **        |                                                                          */
    /* 54.89+                  *    X     *        +                                                                          */
    /*      |                  *          *        |                                                                          */
    /* 54.89+                  ***      ***        +                                                                          */
    /*      |                    ********          |                                                                          */
    /* 54.89+    ******                            +                                                                          */
    /*      |  ***    ***                          |                                                                          */
    /* 54.89+ **  X XX  **                         +                                                                          */
    /*      | *-- 40m ---*                         |                                                                          */
    /* 54.89+ *    X     *                         +                                                                          */
    /*      | ***      ***      ******             |                                                                          */
    /* 54.89+   ********      ***    ***           +                                                                          */
    /*      |                **        **          |                                                                          */
    /* 54.89+                *    X     *          +                                                                          */
    /*      |                *          *          |                                                                          */
    /* 54.89+                ***      ***          +                                                                          */
    /*      |                  ********            |                                                                          */
    /*      ------+--------+--------+--------+-----                                                                           */
    /*          -1.486   -1.485   -1.484   -1.483                                                                             */
    /*                          X                                                                                             */
    /*                                                                                                                        */
    /*  R Log                                                                                                                 */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /* class       : SpatialPointsDataFrame                                                                                   */
    /* features    : 9                                                                                                        */
    /* extent      : -1.485913, -1.482129, 54.89879, 54.90083  (xmin, xmax, ymin, ymax)                                       */
    /* crs         : +proj=longlat +datum=WGS84 +no_defs                                                                      */
    /* variables   : 1                                                                                                        */
    /* names       : ID                                                                                                       */
    /* min values  :  1                                                                                                       */
    /* max values  :  9                                                                                                       */
    /* [1] "mdist"                                                                                                            */
    /*             [,1]      [,2]       [,3]      [,4]       [,5]       [,6]                                                  */
    /*  [1,]   0.000000  11.78930   6.900236  92.63664 281.951880 284.426758                                                  */
    /*  [2,]  11.789302   0.00000  12.175717  82.84245 270.262401 272.735074                                                  */
    /*  [3,]   6.900236  12.17572   0.000000  87.85984 279.560231 282.043526                                                  */
    /*  [4,]  92.636642  82.84245  87.859840   0.00000 201.290456 203.795350                                                  */
    /*  [5,] 281.951880 270.26240 279.560231 201.29046   0.000000   2.506205                                                  */
    /*  [6,] 284.426758 272.73507 282.043526 203.79535   2.506205   0.000000                                                  */
    /*  [7,] 291.953173 280.23529 289.686749 212.17560  11.473582   9.174054                                                  */
    /*  [8,] 264.674749 254.63311 259.880002 172.05249 113.135363 114.390804                                                  */
    /*  [9,] 311.844999 300.32672 308.914026 226.76004  38.771199  37.248840                                                  */
    /*             [,7]     [,8]      [,9]                                                                                    */
    /*  [1,] 291.953173 264.6747 311.84500                                                                                    */
    /*  [2,] 280.235286 254.6331 300.32672                                                                                    */
    /*  [3,] 289.686749 259.8800 308.91403                                                                                    */
    /*  [4,] 212.175602 172.0525 226.76004                                                                                    */
    /*  [5,]  11.473582 113.1354  38.77120                                                                                    */
    /*  [6,]   9.174054 114.3908  37.24884                                                                                    */
    /*  [7,]   0.000000 122.1852  36.85969                                                                                    */
    /*  [8,] 122.185196   0.0000 105.23283                                                                                    */
    /*  [9,]  36.859689 105.2328   0.00000                                                                                    */
    /*                                                                                                                        */
    /* Call:                                                                                                                  */
    /* hclust(d = as.dist(mdist), method = "complete")                                                                        */
    /*                                                                                                                        */
    /* Cluster method   : complete                                                                                            */
    /* Number of objects: 9                                                                                                   */
    /*                                                                                                                        */
    /* [1] "xy$clust"                                                                                                         */
    /* [1] 1 1 1 2 3 3 3 4 3                                                                                                  */
    /*                 min       max                                                                                          */
    /* coords.x1 -1.486913 -1.481129                                                                                          */
    /* coords.x2 54.897790 54.901830                                                                                          */
    /* [1] "cent"                                                                                                             */
    /*           [,1]     [,2]                                                                                                */
    /* [1,] -1.482201 54.90079                                                                                                */
    /* [2,] -1.482880 54.90011                                                                                                */
    /* [3,] -1.485821 54.89927                                                                                                */
    /* [4,] -1.484275 54.89879                                                                                                */
    /* class    : CirclesRange                                                                                                */
    /*                                                                                                                        */
    /* variables: X1 X2                                                                                                       */
    /*                                                                                                                        */
    /*                                                                                                                        */
    /* presence points: 4                                                                                                     */
    /*          X1       X2                                                                                                   */
    /* 1 -1.482201 54.90079                                                                                                   */
    /* 2 -1.482880 54.90011                                                                                                   */
    /* 3 -1.485821 54.89927                                                                                                   */
    /* 4 -1.484275 54.89879                                                                                                   */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

     /*             _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */
