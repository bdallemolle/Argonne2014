// Include standard libraries
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
// Include OpenCV libraries
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv/cv.h>

// Included namspaces
using namespace cv;
using namespace std;

// Globals
bool debug;

/* dub_to_str()
 Converts a double
 to a C++ string */
void dub_to_str(double x, string &str) {
    char c_string[16];
    sprintf(c_string, "%f", x);
    str += c_string;
    return;
}

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

void parse(int argc, char** argv,         // Commandline arguments
           int &lower, int &upper,        // Threshold bounds
           int &inc,                      // Increment value (delta)
           int &xdelta, int &ydelta)
{
    vector< vector<string> > commands;          // List of commands/options
    cmdlineParse(argc, argv, commands);         // Parse arguements
    
    for (int i = 1; i < commands.size(); i++) {
        if (!(strcmp("-thresh", commands[i][0].c_str())) ||
            !(strcmp("-t", commands[i][0].c_str())))  {
            // Read threshold
            if (commands[i].size() < 3) {
                cerr << "Must provide 2 arguements for the threshold bounds..." << endl;
                cerr << "Using defaults..." << endl;
            }
            else {
                lower = atoi(commands[i][1].c_str());
                upper = atoi(commands[i][2].c_str());
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
        else if (!strcmp("-x", commands[i][0].c_str())) {
            // Read x bounds
            if (commands[i].size() < 2) {
                cerr << "Must provide an arguement for the horizontal quantization..." << endl;
                cerr << "Using default..." << endl;
            }
            else {
                xdelta = atoi(commands[i][1].c_str());
            }
        }
        else if (!strcmp("-y", commands[i][0].c_str())) {
            // Read y bounds
            if (commands[i].size() < 2) {
                cerr << "Must provide an arguement for the vertical quantization..." << endl;
                cerr << "Using default..." << endl;
            }
            else {
                ydelta = atoi(commands[i][1].c_str());
            }
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

void get_rect(string line, int &thresh, int &width, int &height,
              double &area, int &x, int &y, double &rot) {
    string token;
    istringstream iss(line);
    int i = 0;
    // Get threshold
    while (getline(iss, token, ' ')) {
        switch (i) {
            case 0:
                thresh = atoi(token.c_str());
                cout << "THRESH: " << thresh << " - ";
                break;
            case 1:
                width = atoi(token.c_str());
                cout << "WIDTH: " << width << " - ";
                break;
            case 2:
                height = atoi(token.c_str());
                cout << "HEIGHT: " << height << " - ";
                break;
            case 3:
                area = atof(token.c_str());
                cout << "AREA: " << area << " - ";
                break;
            case 4:
                x = atoi(token.c_str());
                cout << "X: " << x << " - ";
                break;
            case 5:
                y = atoi(token.c_str());
                cout << "Y: " << y << " - ";
                break;
            case 6:
                rot = atof(token.c_str());
                cout << "ROT: " << rot << endl;
                break;
            default:
                cerr << "ERROR" << endl;
                break;
        }
        i++;
    }
    
    return;
}

// * -------------------------------------------------------------------------------- * //

void bin_results(int** bins, int xdelta, int ydelta,
                 int bin_rows, int bin_cols, Point p) {
    //
    for (int i = 1; i <= bin_rows; i++) {
        if (p.y >= (ydelta * (i - 1)) && p.y < (ydelta * i)) {
            for (int j = 1; j <= bin_cols; j++) {
                if (p.x >= (xdelta * (j - 1))
                    && p.x < (xdelta * j)) {
                    bins[i][j]++;
                }
            }
        }
    }
    return;
};

// * -------------------------------------------------------------------------------- * //

void draw_results(Mat &src, string input_fn, int lower, int upper, int inc,
                  int xdelta, int ydelta, int image_width, int image_height) {
    ifstream inputfile(input_fn);
    string line;
    int thresh = -1;
    int width = -1, height = -1;
    double area = -1;
    int x = -1, y = -1;
    double rot = -1.0;
    int bin_rows = ceil((double)image_height / (double)ydelta);
    int bin_cols = ceil((double)image_width / (double)xdelta);
    int** bins = (int**)malloc(bin_rows * sizeof(int*));
    for (int i = 0; i < bin_rows; i++) bins[i] = (int*)malloc(bin_cols * sizeof(int));
    int count = 0;
    
    cerr << "Bin rows: " << bin_rows << endl;
    cerr << "Bin cols: " << bin_cols << endl;
    
    // Zero out bins
    for (int i = 0; i < bin_rows; i++) {
        for (int j = 0; j < bin_cols; j++) {
            bins[i][j] = 0;
        }
    }
    
    //
    if (inputfile.is_open()) {
        getline(inputfile, line);  // Get first line
        while (getline(inputfile, line)) {
            //
            get_rect(line, thresh, width, height, area, x, y, rot);
            //
            Point p1 = Point(x - floor((0.5)*(double)width), y + floor((0.5)*(double)height));
            Point p2 = Point(x + floor((0.5)*(double)width), y - floor((0.5)*(double)height));
            //
            rectangle(src, p1, p2, Scalar(0, 255, 255), 3);
            // Hist the centers of the rectangles
            bin_results(bins, xdelta, ydelta, bin_rows, bin_cols, Point(x, y));
            count++;
        }
    }
    else {
        cerr << "Unable to open data file, exiting program..." << endl;
        exit(-1);
    }
    
    // Draw numbers on bins
    for (int i = 0; i < bin_rows; i++) {
        for (int j = 0; j < bin_cols; j++) {
            if (bins[i][j] > 0) {
                string str;
                int_to_str(bins[i][j], str);
                // dub_to_str(bins[i][j] / (double)count, str);
                putText(src, str, Point(j * xdelta, i * ydelta), FONT_HERSHEY_PLAIN, 2, Scalar(0, 255, 0), 5);
                // cerr << "Writing text?" << endl;
            }
        }
    }
    
    // Free histogram
    for (int i = 0; i < bin_rows; i++) free(bins[i]);
    free(bins);
    
    return;
}

// * -------------------------------------------------------------------------------- * //

int main(int argc, char** argv) {
    Mat src;                                    // Source image
    int image_width = 0;                        // Image width
    int image_height = 0;                       // Image height
    int xdelta = 150;                            // Delta for histogram x-direction quantization
    int ydelta = 150;                            // Delta for histogram y-direction quantization
    int lower = 5, upper = 255;                 // Filter thresholds
    int inc = 5;                                // Threshold increment bound
    vector< vector<string> > commands;          // List of commands/options
    string image_fn;                        // Input image
    string input_fn;                        // Input filename
    string output_fn;                       // Output filename
    debug = false;                              // Initialize debug output toggle
    
    // Parse input
    parse(argc, argv, lower, upper, inc, xdelta, ydelta);
    
    // Verify filenames exist
    if (argc < 4) {
        cerr << "Must provide input filename and destination filename, exiting..." << endl;
        exit(-1);
    }
    
    // Load source image
    image_fn += argv[1];
    input_fn += argv[2];
    output_fn += argv[3];
    src = imread(image_fn, 1);
    if (!src.data) {
        cerr << "No image data, exiting..." << endl;
        exit(-1);
    }
    
    // Get dimensions of image
    image_width = src.cols;
    image_height = src.rows;
    
    if (debug) {
        cout << "IMAGE: " << image_fn << endl;
        cout << "INPUT DATA: " << input_fn << endl;
        cout << "OUTPUT DATA: " << output_fn << endl;
        
    }
    
    // Output information
    if (debug) {
        cout << "Lower: " << lower << endl;
        cout << "Upper: " << upper << endl;
        cout << "Increment: " << inc << endl;
        cout << endl;
    }
    
    draw_results(src, input_fn, lower, upper, inc,
                 xdelta, ydelta, image_width, image_height);
    
    imwrite(output_fn, src);
    
    return 0;
}