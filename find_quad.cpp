// OpenCV Libraries
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv/cv.h>
// Standard Libraries
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace cv;     // Use OpenCV classes
using namespace std;    // Use standard namespace

bool write_out;         // Write output data to file
bool write_all;         // Write output data and images out
bool over_lay;          // Overlay output over original image
bool debug;             // Print debugging information
bool canny;             // Use Canny edge detection instead of binary threshold

#define N 4             // Define the number of sides on the shape as 4 (quadrilateral)
#define PI 3.14159265   // PI


// UPDATE THIS //
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
/* -angle <int:threshold_increment>                      // Default: 0.5                */
/*   MAXIMUM ANGLE TOLERANCE FROM 90˚ FOR QUAD CORNERS                                  */
/* -x <int:xMin> <int:xMax>                              // Defaults: {0, image.width}  */
/*   LOWER AND UPPER BOUNDS ON RECTANGLE X POSITION                                     */
/* -y <int:yMin> <int:yMax>                              // Defaults: {0, image.height} */
/*   LOWER AND UPPER BOUNDS ON RECTANGLE Y POSITION                                     */
/* -canny <int:lower_thresh>                             // Defaults: 35                */
/*   FLAG TO USE CANNY EDGE DETECTOR AND BOUNDS ON LOWER EDGE THRESHOLD FOR CANNY       */
/* -w                                                                                   */
/*   FLAG TO WRITE RAW OUTPUT TO FILE                                                   */
/* -all                                                                                 */
/*   FLAG TO WRITE RAW OUTPUT AND IMAGES TO FILE                                        */
/* -o                                                                                   */
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

/* find_angle()
 Finds the angle (in degrees) between two points 
 using the Law of Cosines */

double find_angle( Point p1, Point p2, Point p0) {
    // Find distances between points
    double a = sqrt(pow((p0.x - p1.x), 2) + pow((p0.y - p1.y), 2));
    double b = sqrt(pow((p0.x - p2.x), 2) + pow((p0.y - p2.y), 2));
    double c = sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2));
    
    // Compute value of cos(C)
    double cosineC = (c*c - a*a - b*b)/(2*a*b);
    double C = acos(cosineC);
    
    // Convert to degrees
    C = (180 / PI) * C;
    return C;
}

// * -------------------------------------------------------------------------------- * //

/* sortPts()
 Selection sorts the points at the corners of a rectangle (vector size 4)
 Returns clockwise sorted array of points
 (top left -> top right -> bottom right -> bottom left)
 */

void sortPts(vector<Point> &r) {
    int i, j, min;
    
    // Select sort by x coordinates
    for (i = 0; i < N - 1; i++) {
        min = i;
        for (j = i + 1; j < N; j++) {
            if (r[j].x < r[min].x) {
                min = j;
            }
        }
        // Swap min
        swap(r[i], r[min]);
    }
    
    //cout << "SORTED BY X " << r << endl;
    
    // Sort y coordinates relative to x
    if (r[0].y > r[1].y)
        swap(r[0], r[1]);
    if (r[2].y < r[3].y)
        swap(r[2], r[3]);
    
    //cout << "WITH Y? " << r << endl;
    
    return;
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

void parse_input(int argc, char** argv,                 // Commandline arguments
                 int &min_thresh, int &max_thresh,      // Threshold bounds
                 int &inc,                              // Threshold increment
                 double &min_prop, double &max_prop,    // Area (proportionate to image)
                 double &max_angle,                     // Maximum cosine value for quadrilateral angles
                 int &xMin, int &xMax,                  // x boundaries
                 int &yMin, int &yMax,                  // y boundaries
                 int &channel,                          // Color channel to use for edge detection
                 int &canny_low)                        // Lower boundary for Canny edge threshold
{
    vector< vector<string> > commands;          // List of commands/options
    cmdlineParse(argc, argv, commands);         // Parse arguements into command vector
    
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
        else if (!(strcmp("-angle", commands[i][0].c_str())) ||
                 !(strcmp("-e", commands[i][0].c_str()))) {
            // Read max angle
            if (commands[i].size() < 2) {
                cerr << "Must provide an arguements for the +/- angle epsilon..." << endl;
                cerr << "Using default..." << endl;
            }
            else {
                max_angle = atof(commands[i][1].c_str());
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
        else if (!(strcmp("-channel", commands[i][0].c_str())) ||
                 !(strcmp("-ch", commands[i][0].c_str())))   {
            // Read canny threshold
            if (commands[i].size() < 2) {
                cerr << "No channel number provided, using default" << endl;
            }
            if (atoi(commands[i][1].c_str()) < 0 || atoi(commands[i][1].c_str()) > 2) {
                cerr << "Not a valid channel (must be 0, 1, 2)" << endl;
                cerr << "Using default channel" << endl;
            }
            else {
                channel = atoi(commands[i][1].c_str());
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

void find_rect(Mat src_gray,          // Source image
               Mat src,               // Actual source image
               int thresh,            // Threshold parameter
               string output_fn,      // Output filename (for images)
               string rawout_fn,      // Raw data filename
               int min_area,          // Rectangle area lower constraint
               int max_area,          // Rectangle area upper constraint
               double max_angle,      // Maximum cosine value for quadrilateral corners
               int canny_low,         // Lower threshold for Canny edge call
               int xMin, int xMax,    // X-axis location restraints
               int yMin, int yMax)    // Y-axis location restraints
{
    Mat threshold_output;                   // Filtered image destination
    Mat drawing;                            // Rectangle drawing image (for output image)
    Mat overlay;                            // Copy of source image for overlay
    string thresh_str;                      // Threshold (as a string)
    int_to_str(thresh, thresh_str);         // ...convert threshold
    vector<vector<Point> > contours;        // Vector of contour lines
    vector<Vec4i> hierarchy;                // Vector of rectangles (heirarchy?)
    double ratio = 0.0;                     // Ratio of width to height
    int area = 0;                           // Temp area of rectangle
    int x = 0, y = 0;                       // Temp (x, y) of rectangle
    int w = 0, h = 0;                       // Temp width, height of rectangle
    ofstream fileOUT(rawout_fn, ios::app);  // Open <rawout_fn> in append mode
    int n = N;                              // Number of points in square
    vector<vector<Point> >squares;          // Vector of points representing squares
    double cosine[N];                       // Cosines of angles
    double maxCosine;                       //
    double minCosine;                       //
    
    
    if (canny) {
        // Use Canny edge detection
        Canny(src_gray, threshold_output, canny_low, thresh);
        bitwise_not(threshold_output, threshold_output);
        erode(threshold_output, threshold_output, Mat(), Point(-1, -1), 7, 1, 1);
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
    vector<Point> approx;
    
    // Create canvas for drawing rectangles
    if (write_all) {
        drawing = Mat::zeros(threshold_output.size(), CV_8UC3);
        overlay = src.clone();
    }
    
    // Test each contour
    for (size_t i = 0; i < contours.size(); i++) {
        // Approximate contour with accuracy proportional to the contour perimeter
        approxPolyDP(Mat(contours[i]), approx, arcLength(Mat(contours[i]), true) * 0.02, true);
        // Check that polygon is N-sided and convex
        if (approx.size() == N &&
            isContourConvex(Mat(approx))) {
            // Check that quadrilateral is within area bounds
            if (fabs(contourArea(Mat(approx))) > min_area &&
                fabs(contourArea(Mat(approx))) < max_area) {
                // Calculate cosine value for each angle
                maxCosine = 0;
                minCosine = 1;
                
                // Set points
                Point pnt;
                vector <Point> pnts;
                for (int j = 0; j < N; j++) {
                    pnt.x = approx[j].x;
                    pnt.y = approx[j].y;
                    pnts.push_back(pnt);
                }
                
                // Sort points in clockwise order
                sortPts(pnts);
                
                // Loop through
                for (int j = 0; j < N; j++) {
                    Point p0, p1, p2;
                    p0 = pnts[j];
                    if (j - 1 < 0) p1 = pnts[j - 1 + n];
                    else p1 = pnts[j - 1];
                    if (j + 1 >= n) p2 = pnts[j + 1 - n];
                    else p2 = pnts[j + 1];
                    cosine[j] = find_angle(p1, p2, p0);
                    if (debug) cout << "Cosine[" << j << "] " << cosine[j] << endl;
                    maxCosine = max(maxCosine, cosine[j]);
                    minCosine = min(minCosine, cosine[j]);
                    double C = find_angle(p1, p2, p0);
                    
                    // cout << "Law of Cosines method found: " << C << endl;
                }
                if (debug) {
                    cout << "Max Cosine: " << maxCosine;
                    cout << " Min Cosine: " << minCosine << endl;
                    cout << "Cosine Diff: " << maxCosine - minCosine << endl;
                }

                // Check that angles are within error range
                if (maxCosine < max_angle) {
                    squares.push_back(pnts);
                    Scalar color = Scalar(255, 0, 255);       // White lines
                    drawContours(drawing, contours, i, color, 8, 8, vector<Vec4i>(), 0, Point());
                }
            }
        }
    }
    
    size_t num_squares = squares.size();
    for (size_t i = 0; i < num_squares; i++) {
        if (debug) {
            // Print user friendly output
            cout << "THRESH=" << thresh;
            cout << " x0=" << squares[i][0].x << " y0=" << squares[i][0].y;
            cout << " x1=" << squares[i][1].x << " y1=" << squares[i][1].y;
            cout << " x2=" << squares[i][2].x << " y2=" << squares[i][2].y;
            cout << " x3=" << squares[i][3].x << " y3=" << squares[i][3].y;
            cout << endl;
        }
        if (write_out) {
            // Add to raw out file
            // fileOUT << input_fn << " ";
            fileOUT << thresh << " ";
            fileOUT << squares[i][0].x << " " << squares[i][0].y  << " ";
            fileOUT << squares[i][1].x << " " << squares[i][1].y  << " ";
            fileOUT << squares[i][2].x << " " << squares[i][2].y  << " ";
            fileOUT << squares[i][3].x << " " << squares[i][3].y  << " ";
            fileOUT << endl;
        }
        else {
            // Print raw output
            cout << thresh << " ";
            cout << w << " " << h << " ";
            cout << area << " ";
            cout << x << " " << y << endl;
        }
        
        // Write quadrilaterals to image (if outputting images)
        if (write_all) {
            if (!over_lay) {
                polylines(drawing, squares, true, Scalar(0, 255, 0), 3, LINE_AA);
            }
            // Draw rectangle on source image
            else polylines(overlay, squares, true, Scalar(0, 255, 0), 3, LINE_AA);
        }
    }
    
    // Write image (drawing) which contains found rectanges
    if (write_all && num_squares > 0) {
        if (debug) cout << "WRITING IMAGE " << output_fn << endl;
        output_fn += "_" + thresh_str;
        output_fn += ".jpg";
        if (!over_lay) {
            // Drawer contours on drawing
            for (size_t i = 0; i < contours.size(); i++ ) {
                //Scalar color = Scalar(255, 255, 255);       // White lines
                //drawContours(drawing, contours, i, color, 1, 8, vector<Vec4i>(), 0, Point());
            }
            imwrite(output_fn, drawing);
        }
        // Write over the input image
        else imwrite(output_fn, overlay);
    }
    
    fileOUT.close();
    
    return;
}

// * -------------------------------------------------------------------------------- * //

int main( int argc, char** argv ) {
    Mat src;                                // Source image
    Mat src_gray;                           // Source image, grayscale filtered
    int input_height = 0, input_width = 0;  // Image dimensions
    int input_area = 0;                     // Image area
    int min_thresh = 5, max_thresh = 255;   // Filter thresholds
    int inc = 5;                            // Threshold increment bound
    int canny_low = 35;                     // Lower threshold for Canny edge call
    int min_area = 0, max_area = 0;         // Rectangle size constraints
    double min_prop = 0.0, max_prop = 1.0;  // Rectangle area bounds relative to image
    int xMin = 0, xMax = 0;                 // x boundaries
    int yMin = 0, yMax = 0;                 // y boundaries
    double max_angle = 95.0;                // Maximum cosine value for quadrilateral angles
    vector< vector<string> > commands;     // List of commands/options
    string input_fn;     // Input filename (includes directory)
    string output_fn;        // Output filename (includes directory)
    string rawout_fn;        // Output filename (includes directory)
    write_out = false;                     // Initialize write out toggle (raw data only)
    write_all = false;                     // Initialize write all toggle (images & raw data)
    debug = false;                         // Initialize output dump toggle
    canny = false;                         // Initialize canny toggle
    int channel = -1;                      // Color channel to use (default = -1 = grayscale)
    
    parse_input(argc, argv, min_thresh, max_thresh, inc,
                min_prop, max_prop, max_angle, xMin, xMax,
                yMin, yMax, channel, canny_low);
    
    // * ---------- INPUT STUFF -------------- * //
    
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
    rawout_fn += ".txt";
    
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
        cout << "Area relative to image bounds: " << " {" << min_prop << ", " << max_prop << "}" << endl;
        cout << "Area bounds: " << " {" << min_area << ", " << max_area << "}" << endl;
        cout << "Max angle cosine: " << max_angle << " Degrees" << endl;
        cout << "x-axis bounds: {" << xMin << ", " << xMax << "}" << endl;
        cout << "y-axis bounds: {" << yMin << ", " << yMax << "}" << endl;
        cout << endl;
    }
    
    // * ---------- MAIN PROCEDURE -------------- * //
    
    // Convert image to color scale
    if (channel != -1) {
        Mat channels[3];
        split(src, channels);
        src_gray = channels[channel];
    }
    else {
        cvtColor(src, src_gray, CV_BGR2GRAY);
    }
    // Blur image
    // blur(src_gray, src_gray, Size(3,3));
    GaussianBlur(src_gray, src_gray, Size(7, 7), 5);
    
    // Add column headers
    if (write_out) {
        ofstream fileOUT(rawout_fn, ios::trunc);  // Open <rawout_fn> in append mode
        fileOUT << "Threshold x0 y0 x1 y1 x2 y2 x3 y3 " << endl;
        fileOUT.close();
    }
    
    // Loop through thresholds
    for (int i = min_thresh; i < max_thresh; i += inc) {
        find_rect(src_gray, src, i, output_fn, rawout_fn, min_area, max_area,
                  max_angle, canny_low, xMin, xMax, yMin, yMax);
    }
    
    cout << endl;

    return 0;
}

