// OpenCV Libraries
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv/cv.h>
// Standard C++ Libraries
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

using namespace cv;
using namespace std;

// Range structure
struct RANGE {
    int xmin;   // Lower bound of x-bin
    int xmax;   // Upper bound of x-bin
    int ymin;   // Lower bound of y-bin
    int ymax;   // Upper bound of y-bin
    int xavg;   // Average value of x-bin
    int yavg;   // Average value of y-bin
    int xmod;   // Mode of x-bin
    int ymod;   // Mode of y-bin
    int xmed;   // Median of x-bin
    int ymed;   // Median of y-bin
};

// *---------------------------------------------------------------- * //

void sort_points(vector <int> &points) {
    int i, j, min;
    int n = points.size();
    
    // Select sort by x coordinates
    for (i = 0; i < n - 1; i++) {
        min = i;
        for (j = i + 1; j < n; j++) {
            if (points[j] < points[min]) {
                min = j;
            }
        }
        // Swap min
        swap(points[i], points[min]);
    }
    
    return;
}

// *---------------------------------------------------------------- * //

void sort_corners(RANGE r[4]) {
    int i, j, min;
    int n = 4;
    
    // Select sort by x coordinates
    for (i = 0; i < n - 1; i++) {
        min = i;
        for (j = i + 1; j < n; j++) {
            if (r[j].xavg < r[min].xavg) {
                min = j;
            }
        }
        // Swap min
        swap(r[i], r[min]);
    }
    
    // Sort y coordinates relative to x
    if (r[0].yavg > r[1].yavg)
        swap(r[0], r[1]);
    if (r[2].yavg < r[3].yavg)
        swap(r[2], r[3]);
    
    return;
}

// *---------------------------------------------------------------- * //

void do_bin(int* hist,                  //
            vector <vector<int> > pnts, //
            unsigned int low,           //
            unsigned int hi,            //
            int delta,                  //
            int point)                  //
{
    // For each point in points, place in bin
    for (int i = 0; i < pnts.size(); i++) {
        int indx = 0;
        // For each bin, see if point belongs to bin
        for (int j = low; j <= hi; j += delta) {
            if (pnts[i][point] > j && pnts[i][point] <= (j + delta)) {
                hist[indx]++;
            }
            indx++;
        }
    }
    return;
}

// *---------------------------------------------------------------- * //

void get_max(int* hist,
             int &max,
             int &max_indx,
             int nBins,
             int exclude)
{
    // Find the max bin
    max = 0;                //
    max_indx = 0;           //
    for (int i = 0; i < nBins; i++) {
        if (hist[i] > max && i != exclude) {
            max_indx = i;
            max = hist[i];
        }
    }
    return;
}

// *---------------------------------------------------------------- * //

int average_point(vector< vector<int> > pnts,
                  int min, int max, int max_val,
                  int pos)
{
    // Find the average within a point range
    int sum = 0;
    for (int i = 0; i < pnts.size(); i++) {
        if (pnts[i][pos] > min &&
            pnts[i][pos] <= max) {
            sum += pnts[i][pos];
        }
    }
    
    return floor((double)sum / (double)max_val);
}

// *---------------------------------------------------------------- * //

int median_point(vector< vector<int> > pnts,
                  int min, int max, int max_val,
                  int pos)
{
    vector <int> bin_pnts;
    // Find the average within a point range
    for (int i = 0; i < pnts.size(); i++) {
        if (pnts[i][pos] > min &&
            pnts[i][pos] <= max) {
            bin_pnts.push_back(pnts[i][pos]);
        }
    }
    
    sort_points(bin_pnts);
    cout << "Printing sorted points: ";
    for (int j = 0; j < bin_pnts.size(); j++)
        cout << bin_pnts[j] << " ";
    cout << endl;
    
    if (bin_pnts.size() % 2 == 0) {
        cout << "EVEN, FIND AVERAGES FOR MEDIAN" << endl;
        int index1 = floor((double)bin_pnts.size()/2.0);
        int index2 = index1 + 1;
        cout << "Values to average (for median): " << index1;
        cout << " " << index2 << endl;
        cout << "Median: " << floor((double)(bin_pnts[index1] + bin_pnts[index2])/2.0) << endl;
        return floor((double)(bin_pnts[index1] + bin_pnts[index2])/2.0);
    }
    else {
        cout << "ODD, PICK MIDDLE" << endl;
        int index = floor((double)bin_pnts.size()/2.0) + 1;
        cout << "Median: " << bin_pnts[index] << endl;
        return bin_pnts[index];
    }
    
    return -1;
}

// *---------------------------------------------------------------- * //

void output(Mat &out,               // Output image
            RANGE r[4],             // Best corner candidates
            RANGE r2[4],            // Runner up corner candidates
            int corner_mode)        // Which corner method to use
{
    // Draw points on image
    for (int i = 0; i < 4; i++) {
        Point pnt;
        if (corner_mode == 0) {
            pnt.x = r[i].xavg;
            pnt.y = r[i].yavg;
        }
        else if (corner_mode == 1) {
            pnt.x = r[i].xmod;
            pnt.y = r[i].ymod;
        }
        else {
            pnt.x = r[i].xmed;
            pnt.y = r[i].ymed;
        }
        circle(out, pnt, 10, Scalar(255, 255, 255), 25);
        // Draw runner up points
        if (corner_mode == 0) {
            pnt.x = r2[i].xavg;
            pnt.y = r2[i].yavg;
        }
        else if (corner_mode == 1) {
            pnt.x = r2[i].xmod;
            pnt.y = r2[i].ymod;
        }
        else {
            pnt.x = r2[i].xmed;
            pnt.y = r2[i].ymed;
        }
        if (pnt.x > 0 && pnt.y > 0 && pnt.x < 4000 && pnt.y < 3000)
            circle(out, pnt, 10,  Scalar(200, 200, 200), 20);
    }
    return;
}

// *---------------------------------------------------------------- * //

void parse_string(string input,             // Input line
                  vector<string> &line)     // Vector of parsed input
{
    int i, j, low, hi;
    int len = input.length();
    char buffer[80];
    
    j = 0, low = 0;
    for (i = 0; i < len + 1; i++) {
        if ((input[i] == ' ') || (input[i] == '\0')) {
            hi = i;
            input.copy(buffer, hi - low, low);
            buffer[hi - low] = '\0';
            line.push_back(string(buffer));
            low = hi + 1;
            j++;
        }
    }
}

// *---------------------------------------------------------------- * //

void parse_input(int argc, char* argv[],            // Commandline args
                 ifstream &in_file,                 // Input file stream
                 string &input_fn,                  // Input file name
                 string &image_fn,                  // Input image name
                 string &output_fn,                 // Output file name
                 string &imageout_fn,               // Image output name
                 unsigned int &rows,                // Number of rows in image
                 unsigned int &cols,                // Number of cols in image
                 int &delta,                        // Size of hist bins
                 vector <vector <int> > &x_pts,     // Vector of corners, x
                 vector <vector <int> > &y_pts,     // Vector of corners, y
                 int &corner_mode,
                 bool &debug)                       // Debug flag
{
    Mat src;                        // Source image matrix
    vector <int> xtemp, ytemp;      // Temporary coordinate placeholder vectors
    vector <string> line;           // Parsed line of input
    string input;                   // Unparsed line of input
    int cnt = 0;                    // Counter variable
    
    // Check for valid number of commandline args
    if (argc < 4) {
        cerr << "Too few commandline args. Exiting..." << endl;
        exit(-1);
    }
    
    // Collect commandline args
    input_fn += argv[1];
    in_file.open(input_fn);         // Get data file
    cout << "file: " << argv[1] << endl;
    delta = atoi(argv[2]);          // Get delta arg
    image_fn += argv[3];
    src = imread(image_fn, 1);
    if (!src.data) {
        cerr << "No image data, exiting..." << endl;
        cerr << endl;
        exit(-1);
    }
    output_fn += argv[4];
    
    
    // Check for corner mode
    if (argc > 5) {
        corner_mode = atoi(argv[5]);
        imageout_fn += argv[5];
        imageout_fn += "-";
    }
    // Check for debug flag
    if (argc > 6 && !(strcmp(argv[6], "-d"))) {
        debug = true;
    }
    
    imageout_fn += argv[3];
    
    // Get image dimensions
    rows = src.rows;
    cols = src.cols;
    cout << "Rows: " << rows << " Cols: " << cols << endl;
    
    // Get first line
    getline(in_file, input);
    
    // Get other lines
    do {
        getline(in_file, input);
        parse_string(input, line);
        
        if (line.size() >= 9) {
            cnt++;
            // Get vectors
            xtemp.push_back(atoi(line[1].c_str()));
            xtemp.push_back(atoi(line[3].c_str()));
            xtemp.push_back(atoi(line[5].c_str()));
            xtemp.push_back(atoi(line[7].c_str()));
            ytemp.push_back(atoi(line[2].c_str()));
            ytemp.push_back(atoi(line[4].c_str()));
            ytemp.push_back(atoi(line[6].c_str()));
            ytemp.push_back(atoi(line[8].c_str()));
            // Add to larger vector
            x_pts.push_back(xtemp);
            y_pts.push_back(ytemp);
        }
        // Clear lines
        line.clear();
        xtemp.clear();
        ytemp.clear();
        
    } while (!in_file.eof());
    // Check for empty file
    if (cnt == 0) {
        cerr << "No points to classify, exiting..." << endl;
        cerr << endl;
        exit(-1);
    }
    
    return;
}

// *---------------------------------------------------------------- * //

void do_hist(vector <vector<int> > x_pts,  // Vector of vectors of x-coordinates
             vector <vector<int> > y_pts,  // Vector of vectors of y-coordinates
             int delta,                    // Bin size
             RANGE r[4],                   // 1st place bin
             RANGE r2[4],                  // 2nd place bin
             unsigned int xbound,          // Upper bound of x-coordinates
             unsigned int ybound)          // Upper bound of y-coordinates
{
    int indx;                   // Histogram index variable
    int num_bins;               // Number of histogram bins
    int i, j, k;                // Index variables
    int* hist;                  // Histogram array
    int max, max2;              // Max bin counts
    int max_indx, max_indx2;    // Max bin index
    double ratio, ratio2;       // Ratio of bins
    

    num_bins = floor(xbound / delta);
    cout << "#bins " << num_bins << " range " << delta << endl;
    hist = (int*)malloc(num_bins* sizeof(int));
    
    // Index, hist, and calculate average
    for (k = 0; k < 4; k++) {
        // Clear the bins
        for (i = 0; i < num_bins; i++) hist[i] = 0;
        // Place coordinates in bins
        do_bin(hist, x_pts, 0, xbound, delta, k);
        
        // Find the max bin
        get_max(hist, max, max_indx, num_bins, -1);
        // Find the 2nd max bin
        get_max(hist, max2, max_indx2, num_bins, max_indx);
        
        // Calculate max bin bounds
        ratio = (double)max / (double)x_pts.size();
        cout << "Best X Ratio:  " << ratio << endl;
        r[k].xmin = max_indx * delta;
        r[k].xmax = (max_indx + 1) * delta;
        // Calculate max2 bin bounds
        ratio2 = (double)max2 / (double)x_pts.size();
        cout << "Runner-up X Ratio: " << ratio2 << endl;
        r2[k].xmin = max_indx2 * delta;
        r2[k].xmax = (max_indx2 + 1) * delta;
        
        // Find average value of point in best bin
        r[k].xavg = average_point(x_pts, r[k].xmin, r[k].xmax, hist[max_indx], k);
        r[k].xmed = median_point(x_pts, r[k].xmin, r[k].xmax, hist[max_indx], k);
        
        // Find average value of point in 2nd best bin
        r2[k].xavg = average_point(x_pts, r2[k].xmin, r2[k].xmax, hist[max_indx2], k);
    }
    
    for (k = 0; k < 4; k++) {
        // Clear the bins
        for (i = 0; i < num_bins; i++) hist[i] = 0;
        // Place y-coordinates in bins
        do_bin(hist, y_pts, 0, ybound, delta, k);
        
        // Find the max bin
        get_max(hist, max, max_indx, num_bins, -1);
        // Find the 2nd max bin
        get_max(hist, max2, max_indx2, num_bins, max_indx);
        
        // Calculate max bin bounds
        ratio = (double)max / (double)y_pts.size();
        cout << "Best Y Ratio: " << ratio << endl;
        r[k].ymin = max_indx * delta;
        r[k].ymax = (max_indx + 1) * delta;
        // Calculate max2 bin bounds
        ratio2 = (double)max2 / (double)y_pts.size();
        cout << "Runner-up Y Ratio: " << ratio2 << endl;
        r2[k].ymin = max_indx2 * delta;
        r2[k].ymax = (max_indx2 + 1) * delta;
        
        // Find average value of point in best bin
        r[k].yavg = average_point(y_pts, r[k].ymin, r[k].ymax, hist[max_indx], k);
        r[k].ymed = median_point(y_pts, r[k].ymin, r[k].ymax, hist[max_indx], k);
        
        // Find average value of point in 2nd best bin
        r2[k].yavg = average_point(y_pts, r2[k].ymin, r2[k].ymax, hist[max_indx2], k);
        // median_point(y_pts, r2[k].ymin, r2[k].ymax, hist[max_indx], k);
    }

    // Free hist
    free(hist);
    
    hist = (int*)malloc(delta*sizeof(int));
    
    // Find mode of x coordinates
    for (k = 0; k < 4; k++) {
        // Clear the bins
        for (i = 0; i < delta; i++) hist[i] = 0;
        // Place x-coordinates in bins
        do_bin(hist, x_pts, r[k].xmin, r[k].xmax, 1, k);
        // Get best mode
        get_max(hist, max, max_indx, delta, -1);
        r[k].xmod = max_indx + r[k].xmin;
        
        // Clear the bins
        for (i = 0; i < delta; i++) hist[i] = 0;
        // Place x-coordinates in bins
        do_bin(hist, x_pts, r2[k].xmin, r2[k].xmax, 1, k);
        // Get best mode
        get_max(hist, max2, max_indx2, delta, -1);
        r2[k].xmod = max_indx2 + r[k].xmin;
    }
    
    // Find mode of y coordinates
    for (k = 0; k < 4; k++) {
        // Clear the bins
        for (i = 0; i < delta; i++) hist[i] = 0;
        // Place x-coordinates in bins
        do_bin(hist, y_pts, r[k].ymin, r[k].ymax, 1, k);
        // Get best mode
        get_max(hist, max, max_indx, delta, -1);
        r[k].ymod = max_indx + r[k].ymin;
        
        // Clear the bins
        for (i = 0; i < delta; i++) hist[i] = 0;
        // Place x-coordinates in bins
        do_bin(hist, y_pts, r2[k].ymin, r2[k].ymax, 1, k);
        // Get best mode
        get_max(hist, max2, max_indx2, delta, -1);
        r2[k].ymod = max_indx2 + r[k].ymin;
    }
    
    free(hist);
    
    cout << endl;
    
    return;
}

// *---------------------------------------------------------------- * //

int main(int argc, char* argv[]) {
    string input_fn;       // Input data file
    string image_fn;        // Input image
    string output_fn;      // Output data file
    string imageout_fn;    // Output image (for debugging)
    ifstream in_file;                 // Output file stream
    vector <vector <int> > x_pts;     // x-coordinates
    vector <vector <int> > y_pts;     // y-coordinates
    unsigned int rows;             // Number of rows in image
    unsigned int cols;             // Number of cols in image
    int delta;                     // Size of bins
    RANGE r[4];                    // Best corners
    RANGE r2[4];                   // Runner up corners
    int i, j;                      // Index variables
    bool debug = true;            // Debug flag
    int corner_mode = 0;           // 0=Average; 1=Mode; 2=Median
    
    parse_input(argc, argv, in_file, input_fn, image_fn, output_fn,
                imageout_fn, rows, cols, delta, x_pts, y_pts, corner_mode, debug);
    
    // Do binning
    do_hist(x_pts, y_pts, delta, r, r2, cols, rows);
    
    // Sort coordinates in correct order
    sort_corners(r);
    
    for (i = 0; i < 4; i++) {
        cout << "X = Range: {" << r[i].xmin << ", ";
        cout << r[i].xmax << "} ";
        cout << "Avg: " << r[i].xavg << endl;
        cout << "Y = Range: {" << r[i].ymin << ", ";
        cout << r[i].ymax << "} ";
        cout << "Avg: " << r[i].yavg << endl;
    }
    
    for (i = 0; i < 4; i++) {
        cout << i << ": ";
        cout << "(" << r[i].xavg;
        cout << ", " << r[i].yavg << ")" << endl;
        cout << "RUNNER UP: " << r2[i].xavg;
        cout << "(" << r2[i].xavg;
        cout << ", " << r2[i].yavg << ")" << endl;
        cout << "MODE: ";
        cout << "(" << r[i].xmod;
        cout << ", " << r[i].ymod << ")" << endl;
        cout << "RUNNER UP: ";
        cout << "(" << r2[i].xmod;
        cout << ", " << r2[i].ymod << ")" << endl;
    }
    
    // Add column headers
    ofstream fileOUT(output_fn, ios::app);  // Open <rawout_fn> in append mode
    
    for (i = 0; i < 4; i++) {
        fileOUT << r[i].xavg << " ";
        fileOUT << r[i].yavg << " ";
    }
    fileOUT << endl;
    fileOUT.close();
    
    if (debug) {
        Mat out = imread(image_fn, 1);
        output(out, r, r2, corner_mode);
        imwrite(imageout_fn, out);
    }
    
    cout << endl;
    
    return 0;
}
