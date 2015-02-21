// OpenCV Libraries
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv/cv.h>
// Standard Libraries
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace cv;     // Use OpenCV classes
using namespace std;    // Use standard namespace

bool write_out;         // Write output data to file
bool write_all;         // Write output data and images out
bool over_lay;          // Overlay output over original image
bool debug;             // Print debugging information
bool canny;             // Use Canny edge detection instead of binary threshold
bool mass;              // Toggle center of mass error checking

#define PI 3.14159265   // PI

/****************************************************************************************/
/****************************************************************************************/
/* USAGE:                                                                               */
/* ./find_rect <input_filename.jpg> <output_filename> + OPTIONAL FLAGS                  */
/* OPTIONAL FLAGS:                                                                      */
/* -thresh <int:min_thresh> <int:max_thresh>                                            */
/*   LOWER AND UPPER BOUNDS FOR BINARY THRESHOLDING, UPPER THRESHOLD OF CANNY           */
/* -inc <int:threshold_increment>                        // Default: 5                  */
/*   INCREMENT OR BOUND OF THRESHOLD INCREMENT IN MAIN LOOP                             */
/* -area <float:min_prop> <float:max_prop>               // Defaults: {0.0, 1.0}        */
/*   LOWER AND UPPER BOUNDS FOR RECTANGLE AREA PROPORTIONATE TO IMAGE                   */
/* -ratio <float:min_ratio> <float:max_ratio>            // Defaults: {0.0, 2.0}        */
/*   LOWER AND UPPER BOUNDS FOR WIDTH/HEIGH RATIO OF RECTANGLES                         */
/* -angle <int:epsilon>                                  // Default: 5                  */
/*   ROTATION TOLERANCE AS MEASURED BY ANGLE VARIANCE (IN DEGREES)                      */
/* -x <int:xMin> <int:xMax>                              // Defaults: {0, image.width}  */
/*   LOWER AND UPPER BOUNDS ON RECTANGLE X POSITION                                     */
/* -y <int:yMin> <int:yMax>                              // Defaults: {0, image.height} */
/*   LOWER AND UPPER BOUNDS ON RECTANGLE Y POSITION                                     */
/* -canny <int:lower_thresh>                             // Defaults: 35                */
/*   FLAG TO USE CANNY EDGE DETECTOR AND BOUNDS ON LOWER EDGE THRESHOLD FOR CANNY       */
/* -write                                                                               */
/*   FLAG TO WRITE RAW OUTPUT TO FILE                                                   */
/* -all                                                                                 */
/*   FLAG TO WRITE RAW OUTPUT AND IMAGES TO FILE                                        */
/* -overlay                                                                             */
/*   IF WRITING IMAGE OUTPUT, OVERLAY BOXES OVER ORIGINAL IMAGE                         */
/* -debug                                                                               */
/*   TURNS ON DEBUGGING MESSAGES, MORE USER FRIENDLY OUTPUT                             */
/****************************************************************************************/
/****************************************************************************************/

/* int_to_str()
 Converts an positive integer < 4-digits
 or a negative integer < 3-digits
 to a C++ string */
void int_to_str(int x, string &str) {
    char c_string[4];
    sprintf(c_string, "%d", x);
    str += c_string;
    return;
}

// * -------------------------------------------------------------------------------- * //

/* get_distance()
 
 */
double get_distance(Point2f p0, Point2f p1) {
    double x0 = p0.x;
    double y0 = p0.y;
    double x1 = p1.x;
    double y1 = p1.y;
    double dist = sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2));
    return dist;
}

// * -------------------------------------------------------------------------------- * //

/* get_vectortip()
 
 */
Point2f get_vectortip(Point2f p0,
                   double theta_degrees,
                   double r)
{
    Point2f p1;
    double theta_radians = theta_degrees * (PI / 180);
    double a = r * cos(theta_radians);
    double b = r * sin(theta_radians);
    p1.x = p0.x + a;
    p1.y = p0.y + b;
    return p1;
}

// * -------------------------------------------------------------------------------- * //

/* cmdlineParse()
 Commandline arguement parser for C++
 cmdlineParse() -> returns vector of vectors of strings by reference
 - Parses arguements based with dash (-) flags
 */
void cmdlineParse(int argc, char** argv,
                  vector< vector<string> > &commands){
    int arg_cnt = 0;
    vector<string> temp;
    
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            // New flag
            temp.clear();
            temp.push_back(argv[i]);
            // Add flag
            commands.push_back(temp);
        }
        else {
            // Arguement
            if (i == 1) {
                // First arguement case (no previous flag)
                temp.clear();
                temp.push_back(argv[i]);
                commands.push_back(temp);
            }
            else commands.back().push_back(argv[i]);
        }
    }
    return;
}

// * -------------------------------------------------------------------------------- * //

/* sortPts()
 Selection sorts the points at the corners of a rectangle (array size 4)
 Returns clockwise sorted array of points 
 (top left -> bottom left -> bottom right -> top left)
 */
void sortPts(Point2f r[4]) {
    int i, j, min;
    int n = 4;
    
    // Select sort by x coordinates
    for (i = 0; i < n - 1; i++) {
        min = i;
        for (j = i + 1; j < n; j++) {
            if (r[j].x < r[min].x) {
                min = j;
            }
        }
        // Swap min
        swap(r[i], r[min]);
    }
    
    // Sort y coordinates relative to x
    if (r[0].y > r[1].y)
        swap(r[0], r[1]);
    if (r[2].y < r[3].y)
        swap(r[2], r[3]);
    
    return;
}

// * -------------------------------------------------------------------------------- * //

void parse_input(int argc, char** argv,                 // Commandline arguments
                 int &min_thresh, int &max_thresh,      // Threshold bounds
                 int &inc,                              // Increment value (delta)
                 double &min_prop, double &max_prop,    // Area (relative to image) bounds
                 double &min_ratio, double &max_ratio,  // Aspect ratio bounds
                 int &epsilon,                          // Angle varience, rotation tolerance (epsilon)
                 double &radius,                           //
                 int &xMin, int &xMax,                  // x-coordinate bounds
                 int &yMin, int &yMax,                  // y-coordinate bounds
                 int &canny_low)                        // Canny threshold lower bound
{
    vector< vector<string> > commands;          // List of commands/options
    cmdlineParse(argc, argv, commands);
    
    for (int i = 1; i < commands.size(); i++) {
        if (!(strcmp("-thresh", commands[i][0].c_str())) ||
            !(strcmp("-t", commands[i][0].c_str())))  {
            // Read threshold
            if (commands[i].size() < 3) {
                cerr << "Must provide 2 arguements for the threshold bounds..." << endl;
                cerr << "Using defaults..." << endl;
            }
            else {
                min_thresh = atoi(commands[i][1].c_str());
                max_thresh = atoi(commands[i][2].c_str());
            }
        }
        else if (!(strcmp("-inc", commands[i][0].c_str())) ||
                 !(strcmp("-i", commands[i][0].c_str())))   {
            // Read increment
            if (commands[i].size() < 2) {
                cerr << "Must provide an arguements for the increment..." << endl;
                cerr << "Using default..." << endl;
            }
            else {
                inc = atoi(commands[i][1].c_str());
            }
        }
        else if (!(strcmp("-area", commands[i][0].c_str())) ||
                 !(strcmp("-a", commands[i][0].c_str()))) {
            // Read area
            if (commands[i].size() < 3) {
                cerr << "Must provide 2 arguements for the area bounds..." << endl;
                cerr << "Using defaults..." << endl;
            }
            else {
                min_prop = atof(commands[i][1].c_str());
                max_prop = atof(commands[i][2].c_str());
            }
        }
        else if (!(strcmp("-ratio", commands[i][0].c_str())) ||
                 !(strcmp("-r", commands[i][0].c_str()))) {
            // Read ratio
            if (commands[i].size() < 3) {
                cerr << "Must provide 2 arguements for the ratio bounds..." << endl;
                cerr << "Using defaults..." << endl;
            }
            else {
                min_ratio = atof(commands[i][1].c_str());
                max_ratio = atof(commands[i][2].c_str());
            }
        }
        else if (!(strcmp("-angle", commands[i][0].c_str())) ||
                 !(strcmp("-e", commands[i][0].c_str()))) {
            // Read max angle
            if (commands[i].size() < 2) {
                cerr << "Must provide an arguements for the +/- angle epsilon..." << endl;
                cerr << "Using default..." << endl;
            }
            else {
                epsilon = atoi(commands[i][1].c_str());
            }
        }
        else if (!strcmp("-x", commands[i][0].c_str())) {
            // Read x bounds
            if (commands[i].size() < 3) {
                cerr << "Must provide 2 arguements for the horizontal bounds..." << endl;
                cerr << "Using defaults..." << endl;
            }
            else {
                xMin = atoi(commands[i][1].c_str());
                xMax = atoi(commands[i][2].c_str());
            }
        }
        else if (!strcmp("-y", commands[i][0].c_str())) {
            // Read y bounds
            if (commands[i].size() < 3) {
                cerr << "Must provide 2 arguements for the vertical bounds..." << endl;
                cerr << "Using defaults..." << endl;
            }
            else {
                yMin = atoi(commands[i][1].c_str());
                yMax = atoi(commands[i][2].c_str());
            }
        }
        else if (!(strcmp("-canny", commands[i][0].c_str())) ||
                 !(strcmp("-c", commands[i][0].c_str())))   {
            // Read canny threshold
            if (commands[i].size() < 2) {
                cerr << "Using default Canny lower threshold" << endl;
            }
            else {
                canny_low = atoi(commands[i][1].c_str());
            }
            // Canny flag
            canny = true;
        }
        else if (!(strcmp("-mass", commands[i][0].c_str())) ||
                 !(strcmp("-m", commands[i][0].c_str()))) {
            // Read center of mass deviance
            if (commands[i].size() < 2) {
                cerr << "Using default distance radius (VALUE)" << endl;
            }
            else {
                radius = atof(commands[i][1].c_str());
                mass = true;
            }
        }
        else if (!(strcmp("-write", commands[i][0].c_str())) ||
                 !(strcmp("-w", commands[i][0].c_str())))   {
            // Write out flag
            write_out = true;
        }
        else if (!(strcmp("-all", commands[i][0].c_str())))  {
            // Write out flag
            write_out = true;
            // Write images flag
            write_all = true;
        }
        else if (!(strcmp("-overlay", commands[i][0].c_str())) ||
                 !(strcmp("-o", commands[i][0].c_str())))  {
            // Overlay flag
            over_lay = true;
        }
        else if (!(strcmp("-debug", commands[i][0].c_str())) ||
                 !(strcmp("-d", commands[i][0].c_str())))  {
            // Debug flag
            debug = true;
        }
    }
    
    return;
}

// * -------------------------------------------------------------------------------- * //

void rect_found(int i, Point2f rect_points[],
                vector<RotatedRect> minRect,
                vector<vector<Point> > contours,
                Mat drawing, Mat overlay,
                int w, int h, int area, int thresh,
                double radius,
                ofstream &fileOUT) {
    Scalar color;
    int error_dist;
    // Write rectangles to image (if outputting images)
    if (write_all) {
        Point2f direction = get_vectortip(minRect[i].center, minRect[i].angle, 100.0);
        if (!over_lay) {
            // Find center of mass of contour
            Moments mu = moments(contours[i], false);
            Point2f center_of_mass = Point2f(mu.m10 / mu.m00, mu.m01 / mu.m00);
            error_dist = get_distance( center_of_mass, minRect[i].center);
            // Do center of mass error checking
            if (mass) {
                if (error_dist > radius) {
                    cerr << "CENTER OF MASS IS OFF, NO RECT" << endl;
                    cerr << "ERROR DISTANCE = " << error_dist << endl;
                    return;
                }
                else {
                    cerr << "ERROR DISTANCE = " << error_dist << endl;
                    color = Scalar(255, 0, 255);
                    circle( drawing, center_of_mass, 5, color, 2, 8);
                }
            }
            
            /* Find convexity defects */
            // Get convexity hull
            vector<int> hull;
            vector<Point> hull_pnt;
            convexHull(Mat(contours[i]), hull, true);
            convexHull(Mat(contours[i]), hull_pnt, false);
            cout << "Number of points in hull array: " << hull_pnt.size() << endl;
            /**
            for (int j = 0; j < hull_pnt.size(); j++) {
                error_dist = get_distance( hull_pnt[j], minRect[i].center);
                cerr << "ERROR DISTANCE = " << error_dist << endl;
                if (error_dist < 300) {
                    cerr << "CONVEX HULL ERROR" << endl;
                    color = Scalar(0, 0, 255);
                    polylines(drawing, hull_pnt, true, color, 3);
                    color = Scalar(255, 0, 255);
                    circle( drawing, hull_pnt[j], 5, color, 2, 8);
                    return;
                }
                else {
                    // cerr << "HULL ERROR DISTANCE = " << error_dist << " INDEX " << j << endl;
                    color = Scalar(255, 0, 255);
                    polylines(drawing, hull_pnt, true, color, 3);
                }
            }
             **/
            
            // Get convexity defects
            vector<Vec4i> defects;
            convexityDefects(Mat(contours[i]), hull, defects);
            // cerr << "Defects for contour[" << i << "]" << endl;
            for (int j = 0; j < defects.size(); j++) {
                int defectPtIdx = defects[j].val[2];
                // cerr << contours[i][defectPtIdx] << endl;
                color = Scalar(255, 255, 0);
                // circle( drawing, contours[i][defectPtIdx], 2, color, 2, 8);
            }
            
            // Draw rectangle on blank canvas
            color = Scalar(0, 255, 255);        // Yellow Lines
            for (int j = 0; j < 4; j++) {
                line(drawing, rect_points[j], rect_points[(j+1)%4], color, 2, 8 );
            }
            color = Scalar(0, 255, 0);       // Green lines = contour
            drawContours(drawing, contours, i, color, 2, 8, vector<Vec4i>(), 0, Point());
            
            // Draw center of rectangle on image
            color = Scalar(255, 0, 0);
            circle(drawing, minRect[i].center, 5, color, 1);
            color = Scalar(0, 0, 255);
            line(drawing, minRect[i].center, direction, color, 2, 8 );
        }
        else {
            // Draw rectangle on source image
            color = Scalar(0, 255, 255);        // Green Lines
            for (int j = 0; j < 4; j++) {
                line(overlay, rect_points[j], rect_points[(j+1)%4], color, 8, 8 );
            }
            // Draw center of rectangle on image
            color = Scalar(255, 0, 0);
            circle(overlay, minRect[i].center, 5, color, 1);
            color = Scalar(0, 0, 255);
            line(overlay, minRect[i].center, direction, color, 2, 8 );
        }
    }
    
    // Output information
    if (debug) {
        // Print user friendly output
        cout << "THRESH=" << thresh;
        cout << " W=" << w << " H=" << h;
        cout << " ROT: " << fabs(minRect[i].angle);
        cout << " A=" << area << endl;
        cout << "Points at ";
        cout << rect_points[0] << " ";
        cout << rect_points[1] << " ";
        cout << rect_points[2] << " ";
        cout << rect_points[3] << endl;
    }
    if (write_out) {
        // Add to raw out file
        fileOUT << thresh << " ";
        fileOUT << w << " " << h << " ";
        fileOUT << area << " ";
        fileOUT << minRect[i].center.x << " " << minRect[i].center.y;
        fileOUT << " " << fabs(minRect[i].angle) << endl;
    }
    else {
        // Print raw output
        cout << thresh << " ";
        cout << w << " " << h << " ";
        cout << area << " ";
        cout << minRect[i].center.x << " " << minRect[i].center.y << endl;
    }
    
    return;
}

// * -------------------------------------------------------------------------------- * //

void find_rotrect(Mat src_gray,          // Source image
                  Mat src,               // Actual source image
                  int thresh,            // Threshold parameter
                  string output_fn,      // Output filename (for images)
                  string rawout_fn,      // Raw data filename
                  int min_area,          // Rectangle area lower constraint
                  int max_area,          // Rectangle area upper constraint
                  int epsilon,           // Angle varience, rotation tolerance (epsilon)
                  double radius,         //
                  int canny_low,         // Lower threshold for Canny edge call
                  double min_ratio,      // Rectangle w/h ratio lower constraint
                  double max_ratio,      // Rectangle w/h ratio upper constraint
                  int xMin, int xMax,    // X-axis location restraints
                  int yMin, int yMax)    // Y-axis location restraints
{
    Mat threshold_output;                   // Filtered image destination
    Mat drawing;                            // Rectangle drawing image (for output image)
    Mat overlay;                            // Copy of source image for overlay
    string thresh_str;                      // Threshold (as a string)
    int_to_str(thresh, thresh_str);         // ... Convert threshold to string
    ofstream fileOUT(rawout_fn, ios::app);  // Open <rawout_fn> in append mode
    vector<vector<Point> > contours;        // Vector of contour lines
    vector<Vec4i> hierarchy;                // Vector of rectangles (heirarchy?)
    Point2f rect_points[4];                 // Temp rectangle points
    double ratio = 0.0;                     // Ratio of width to height
    int area = 0;                           // Temp area of rectangle
    int x = 0, y = 0;                       // Temp (x, y) of rectangle
    int w = 0, h = 0;                       // Temp width, height of rectangle
    double abs_angle = 0.0;                 // Absolute value of rotation angle
    int cnt = 0;                            // Counter (for rectangles passing filter)
    int i, j;                               // Index variables
    int c1, c2, c3;
    
    if (canny) {
        // Use Canny edge detection
        Canny(src_gray, threshold_output, canny_low, thresh);
        //bitwise_not(threshold_output, threshold_output);
        //erode(threshold_output, threshold_output, Mat(), Point(-1, -1), 1, 1, 1);
    }
    else {
        // Detect edges using BINARY Threshold
        threshold(src_gray, threshold_output, thresh, 255, THRESH_BINARY);
    }
    
    // Find contours
    findContours(threshold_output, contours, hierarchy, RETR_TREE,
                 CHAIN_APPROX_SIMPLE, Point(0, 0));
    
    // Approximate contours to polygons + get bounding rects
    vector< vector<Point> > contours_poly(contours.size());
    vector<RotatedRect> minRect(contours.size());
    
    //
    for (int i = 0; i < contours.size(); i++) {
        minRect[i] = minAreaRect(Mat(contours[i]));
    }
    
    // Create canvas for drawing rectangles
    if (write_all) {
        drawing = Mat::zeros(threshold_output.size(), CV_8UC3);
        overlay = src.clone();
    }
    
    for (i = 0; i < contours.size(); i++) {
        // Get rectangle corner points
        minRect[i].points(rect_points);
        
        // Need to sort points in these rectangles
        sortPts(rect_points);
        x = floor(rect_points[0].x);
        y = floor(rect_points[0].y);
        
        // Get width, height
        h = get_distance(rect_points[0], rect_points[1]);
        w = get_distance(rect_points[0], rect_points[3]);
        // Get ratio
        ratio = (double)w/(double)h;
        // Get area
        area = w * h;
        
        // *** Filter rectangles *** //
        // Check boundary restraints
        if ((x > xMin) && (x < xMax) && (y > yMin) && (y < yMax) &&
            ((x + w) < xMax) && ((y + h) < yMax)) {
            // Check area restraints
            if ((area > min_area) && (area < max_area)) {
                // Check for ratio restraints
                if ((ratio > min_ratio) && (ratio < max_ratio)) {
                    abs_angle = fabs(minRect[i].angle);
                    // Check for rotation angle
                    if ((abs_angle < epsilon) || ((abs_angle < (90 + epsilon))
                         && (abs_angle > (90 - epsilon)))) {
                        // *** Rectangle passed filter! *** //
                        cnt++;  // Increment temp count
                        // Write rectangles to image (if outputting images)
                        rect_found(i, rect_points, minRect, contours,
                                   drawing, overlay, w, h, area, thresh, radius, fileOUT);
                    }
                }
            }
        }
    }
    
    // Write image (drawing) which contains found rectanges
    if (write_all && cnt > 0) {
        if (debug) cout << "WRITING IMAGE " << output_fn << endl;
        if (!over_lay) {
            // Drawer contours on drawing
            for (size_t i = 0; i < contours.size(); i++ ) {
                Scalar color = Scalar(50, 50, 50);       // White lines
                drawContours(drawing, contours, i, color, 1, 8, vector<Vec4i>(), 0, Point());
            }
            imwrite(output_fn, drawing);
        }
        else {
            // Write over the input image
            imwrite(output_fn, overlay);
        }
    }
    
    fileOUT.close();
    
    return;
}

// * -------------------------------------------------------------------------------- * //

int main(int argc, char** argv) {
    Mat src;                                    // Source image
    Mat src_gray;                               // Source image, grayscale filtered
    int input_height = 0, input_width = 0;      // Image dimensions
    int input_area = 0;                         // Image area
    int min_thresh = 5, max_thresh = 255;       // Filter thresholds
    int inc = 5;                                // Threshold increment bound
    int canny_low = 35;                         // Lower threshold for Canny edge call
    int min_area = 0, max_area = 0;             // Rectangle size constraints
    double min_prop = 0.0, max_prop = 1.0;      // Rectangle w/h ratio constraints
    int xMin = 0, xMax = 0;                     // x boundaries
    int yMin = 0, yMax = 0;                     // y boundaries
    double min_ratio = 0.0, max_ratio = 5.0;    // Maximum w/h ratio
    int epsilon = 10;                           // Angle varience, rotation tolerance (epsilon)
    double radius = 100;                           // Radius of center of mass error tolerance
    vector< vector<string> > commands;          // List of commands/options
    string input_fn;                            // Input filename (includes directory)
    string output_fn;                           // Output filename (includes directory)
    string rawout_fn;                           // Output filename (includes directory)
    write_out = false;                          // Initialize write out toggle (raw data only)
    write_all = false;                          // Initialize write all toggle (images & raw data)
    over_lay = false;                           // Initialize overlay output toggle
    debug = false;                              // Initialize debug output toggle
    canny = false;                              // Initialize canny toggle
    mass = false;
    
    parse_input(argc, argv, min_thresh, max_thresh, inc,
                min_prop, max_prop, min_ratio, max_ratio,
                epsilon, radius,
                xMin, xMax, yMin, yMax, canny_low);
    
    // *** INPUT STUFF *** //
    
    // Verify filenames exist
    if (argc < 3) {
        cerr << "Must provide input filename and destination filename, exiting..." << endl;
        exit(-1);
    }
    
    // Load source image
    input_fn += argv[1];
    output_fn += argv[2];
    src = imread(input_fn, 1);
    if (!src.data) {
        cerr << "No image data, exiting..." << endl;
        exit(-1);
    }
    
    if (debug) cout << "OPENING: " << input_fn << endl;
    
    // Establish output filename
    rawout_fn += argv[2];
    
    // Get dimensions and area of image
    input_width = src.cols;
    input_height = src.rows;
    input_area = input_width * input_height;
    
    // Calculate proportionate area
    min_area = (int)(min_prop * input_area);
    max_area = (int)(max_prop * input_area);
    
    // Check for default {x, y} range
    if (xMax == 0) xMax = input_width;
    if (yMax == 0) yMax = input_height;
    
    // Output information
    if (debug) {
        cout << "Image width: " << input_width << endl;
        cout << "Image height: " << input_height << endl;
        cout << "Image area: " << input_area << endl;
        cout << "Threshold: " << " {" << min_thresh << ", " << max_thresh << "}" << endl;
        cout << "Delta/Increment: " << inc << endl;
        cout << "Size proportions: " << " {" << min_prop << ", " << max_prop << "}" << endl;
        cout << "Ratio: " << " {" << min_ratio << ", " << max_ratio << "}" << endl;
        cout << "x-axis bounds: {" << xMin << ", " << xMax << "}" << endl;
        cout << "y-axis bounds: {" << yMin << ", " << yMax << "}" << endl;
        cout << "Angle Epsilon: " << epsilon << endl;
        cout << "Center of Mass Error: " << radius << endl;
        cout << "Area: " << " {" << min_area << ", " << max_area << "}" << endl;
        cout << endl;
    }
    
    // *** MAIN PROCEDURE *** //
    
    // Convert image to gray and blur it
    cvtColor(src, src_gray, CV_BGR2GRAY);
    blur(src_gray, src_gray, Size(3,3));
    
    // Add column headers
    if (write_out) {
        ofstream fileOUT(rawout_fn, ios::trunc);  // Open <rawout_fn> in append mode
        fileOUT << "Threshold Width Height Area x y Rotation" << endl;
        fileOUT.close();
    }
    
    // Loop through thresholds
    for (int i = min_thresh; i < max_thresh; i += inc) {
        if (debug) cout << "Iter: " << i << endl;
        find_rotrect(src_gray, src, i, output_fn, rawout_fn, min_area, max_area,
                  epsilon, radius, canny_low, min_ratio, max_ratio, xMin, xMax, yMin, yMax);
    }
    
    cout << endl;

    return 0;
}

