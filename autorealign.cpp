/* Program to test several preprocessing techniques
 for prepare insect drawer image for rectangle recognition
 By Bryan Dalle Molle
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;
using namespace std;

struct rect {
    Point2f p[4];
    int p_cnt;
};

rect image_rect;
rect ref_rect;

// ------------------------------------------------------------------------------------ //

/* bool parse_cmdline
 Parses commandline, returning appropriate values by reference
 Returns true if valid usage, returns false if invalid
 */

int parse_cmdline(int argc, char** argv,
                  Mat &src, string &input_fn,
                  double xdrawer, double ydrawer,
                  double r, string &output_fn,
                  bool &write)
{
    int i = 0;          // Index variables
    int x[4];           // x-coordinates of corners
    int y[4];           // y-coordinates of corners
    int cols, rows;     // Rows and columns of image
    double k;           // Scaling factor
    int w, h;           // Width, height of rectified rectangle
    int d0;             // x axis offset
    int d1;             // y axis offset
    
    // Check for min args
    if (argc < 10) {
        fprintf(stderr, "\nToo few args (%d)", argc);
        fprintf(stderr, "\nExiting...\n");
        exit(-1);
    }
    
    // Open source image
    input_fn += argv[1];
    src = imread(input_fn);
    output_fn += argv[2];
    // Check for valid image
    if (!src.data) {
        fprintf(stderr, "\nNo image data");
        fprintf(stderr, "\nExiting...\n");
        exit(-1);
    }
    cout << "output fn " << output_fn << endl;
    
    // Register transform reference points
    for (i = 0; i < 4; i++) {
        x[i] = atoi(argv[(2 * i) + 3]);
        y[i] = atoi(argv[(2 * i) + 4]);
        image_rect.p_cnt++;
    }
    for (i = 0; i < 4; i++) {
        image_rect.p[i].x = x[i];
        image_rect.p[i].y = y[i];
    }
    
    // Grab reference dimensions
    if (argc > 12) {
        xdrawer = atof(argv[11]);
        ydrawer = atof(argv[12]);
    }
    else {
        cerr << "No reference dimensions given! Exiting..." << endl;
        exit(-1);
    }
    
    // Grab write reference rectangle flag
    if (argc > 13 && !strcmp(argv[13], "-w")) {
        write = true;
    }
    
    // Calculate rectified drawer coordinates
    cols = src.cols;
    rows = src.rows;
    // Get h, w, and scaling factor k
    h = floor(r * (double)rows);
    k = (double)h / ydrawer;
    w = floor(k * xdrawer);
    // Debug output
    cout << "h = " << h << endl;
    cout << "w = " << w << endl;
    cout << "k = " << k << endl;
    // Get offsets d0 and d1
    d0 = floor(0.5 * (double)(cols - w));
    cout << "d0 = " << d0 << endl;
    d1 = floor(0.5 * (double)(rows - h));
    cout << "d1 = " << d1 << endl;
    /* Get reference rectangle points */
    ref_rect.p[0].x = d0;       // p0
    ref_rect.p[0].y = d1;
    cout << "p0 = " << ref_rect.p[0] << endl;
    ref_rect.p[1].x = d0;       // p1
    ref_rect.p[1].y = d1 + h;
    cout << "p1 = " << ref_rect.p[1] << endl;
    ref_rect.p[2].x = d0 + w;   // p2
    ref_rect.p[2].y = d1 + h;
    cout << "p2 = " << ref_rect.p[2] << endl;
    ref_rect.p[3].x = d0 + w;   // p3
    ref_rect.p[3].y = d1;
    cout << "p3 = " << ref_rect.p[3] << endl;

    return 0;
}

// ------------------------------------------------------------------------------------ //

int main(int argc, char** argv) {
    // Variable dictionary
    Mat src;                        // Source image
    Mat trans_mat;                  // Transformation matrix
    Mat dst;                        // Destination?
    image_rect.p_cnt = 0;           // Input image rectangular coordinates
    ref_rect.p_cnt = 0;             // Reference image rectangular coordinates
    string input_fn;                // Input image filename
    string output_fn;               // Temp output
    //double xdrawer = 21.4325;     // Aspect ratio of museum drawer...
    //double ydrawer = 17.25;
    double ydrawer = 10.4325;       // Aspect ratio of test drawer...
    double xdrawer = 11.5;
    double r = 0.8;                 // Ratio of rectified drawer height to image height
    int i = 0;                      // Index
    bool write = false;
    
    // Parse commandline
    parse_cmdline(argc, argv, src, input_fn, xdrawer, ydrawer, r, output_fn, write);
    
    // Verify input
    cout << "Drawer coordinates for transform: " << endl;
    for (i = 0; i < 4; i++)
        cout << "(" << i << "): " << image_rect.p[i] << endl;
    
    // Get transformation
    trans_mat = getPerspectiveTransform(image_rect.p, ref_rect.p);
    
    // Do warp
    warpPerspective(src, dst, trans_mat, src.size());
    
    if (write) {
        // Rect out_rect = Rect(ref_rect)
        rectangle(dst, ref_rect.p[0], ref_rect.p[2], Scalar(255, 0, 255), 5);
    }

    // Save output image
    imwrite(output_fn, dst);
    
    // Exit
    cout << "\nProgram executed normally. Exiting...\n" << endl;
    
    return 0;
}
