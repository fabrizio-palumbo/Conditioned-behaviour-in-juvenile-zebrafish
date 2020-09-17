#include <thread>
#include "tracker.h"

using namespace cv;
using namespace std;

int const max_BINARY_value = 255;
bool trackObjects = true;
bool useMorphOps = true;
bool condition=false;// boolean to keep track if you are conditioning the animal or not.
const int MAX_NUM_OBJECTS=5;//max number of objects to be detected in frame
int MIN_OBJECT_AREA = 1;//minimum object area in pixel, you can here control for the area of the object detected.
int MAX_OBJECT_AREA = 655;//maximum object area in pixel,you can here control for the area of the object detected.
int color_tresh_max=180;//max object color in a scale 0-255 (black to white), you can here control for the brightnes of the object detected.
int color_tresh_min=1;//minimum object color in a scale 0-255 (black to white), you can here control for the brightnes of the object detected.
int DISTANCE_LIMIT=300;//distance expressed in pixel, it represent the maximum reasonable distance the animal can move across two consecutive detection point. This i a control to avoid tracking noise happening on the side of the arena.

void tracker::init()
{}
void tracker::start()
{}

void tracker::start_tracking(cv::Mat *resultROI, Mat *threshROI, Mat *frameROI, Mat *backROI, Mat *foreROI, BackgroundSubtractor* bgROI, int *xposROI, int *yposROI, int *frame_processedROI)
{
    bool thresh_active;
        double areabig=0;
        QVector <double> area ;
        QVector<Moments> moment ;
        Moments momentbig;
        emit send_clock(false);
        if(*frame_processedROI ==0 ){
            emit send_position(*xposROI, *yposROI , false);
        }
        bgROI->apply( *frameROI,*foreROI);
        //cv::imshow("fore",*foreROI ); use it to display the image of the foreground detected
        if(*frame_processedROI < 300 ){
            /*it is used as a control on animal max displacement across two consecutive frames.
            if no movement is detected within this threshold (300 frames) then this control is removed until further animal detection*/
            thresh_active=true;
        }else
        {thresh_active=false;}
        /* This can be used to speed up the tracking if you have a very stable background, then there is no need for having a dynamic background calculation
if(*frame_processedROI < 400 && !stopback){
                    bgROI->getBackgroundImage(*backROI);
    }*/
        /*resultROI=*foreROI;// commented to reduce computation time, uncomment if you wanna store all the calculated foreground
        threshold( *resultROI, *threshROI, threshold_value, max_BINARY_value,threshold_type ); //thresholding the image can help in case of very noise situation
        cv::imshow("thresh", *threshROI); //it display the thresholded image*/
        Mat erodeElement = getStructuringElement( MORPH_RECT,Size(3,3));
        Mat dilateElement = getStructuringElement( MORPH_RECT,Size(8,8));
        erode(*foreROI,*foreROI,erodeElement);
        dilate(*foreROI,*foreROI,dilateElement);
        if(trackObjects){
            int xtemp,ytemp;
            vector< vector<Point> > contours;
            vector<Vec4i> hierarchy;
            //find contours of filtered image using openCV findContours function
            findContours(*foreROI,contours,hierarchy,CV_RETR_LIST ,CV_CHAIN_APPROX_SIMPLE );
            //use moments method to find our filtered object
            int objectFound = 0;
            int numObjects = contours.size();
            //if number of objects greater than MAX_NUM_OBJECTS we have a noisy filter
            if(numObjects<MAX_NUM_OBJECTS && numObjects>0){
                //std::cout <<hierarchy[index][0]<<endl;
                // std::cout <<contours.size()<<endl;
                for (int i = 0; i<numObjects ; i++) {
                    moment.push_back(moments((cv::Mat)contours[i]));
                    area.push_back(moment[i].m00);
                    //std::cout << area [i]<<endl; if you want to monitor real time the area of the object detected
                }
                if(numObjects>1){
                    //std::cout <<numObjects<<endl; //to monitor how many object you detect at every iteration
                    for (int i = 0; i<numObjects ; i++) {
                        if(i==0){
                            if(area[0]>MIN_OBJECT_AREA && area[0]<MAX_OBJECT_AREA  ){
                                areabig=area[0];
                                momentbig= moment[0];
                            }
                        }else{
                            if(areabig<=area[i] && area[i]>MIN_OBJECT_AREA && area[i]<MAX_OBJECT_AREA  ){
                                areabig=area[i];
                                momentbig= moment[i];
                            }
                        }
                    }
                }
                else {areabig=area[0];
                    momentbig=moment[0];
                }
                if(areabig>MIN_OBJECT_AREA && areabig<MAX_OBJECT_AREA ){
                    xtemp= momentbig.m10/areabig;
                    ytemp= momentbig.m01/areabig;
                   /* //This piece of code can be used if you want to control also for the color intensity of the object you detect
                    // for instance if you know the animal is always gonna be a dark object in the image you can check that
                    // the detected object is not brighter that 150 (on a scale 0-255)
                    //uchar min_col=( frameROI->at<uchar>(ytemp,xtemp));
                    // uchar max_col=(frameROI->at<uchar>(ytemp,xtemp));
                    //if((int)max_col <color_tresh_max && (int)min_col >color_tresh_min ){*/
                    if(thresh_active){
                        if (abs(xtemp-*xposROI)<DISTANCE_LIMIT && abs(ytemp-*yposROI)<DISTANCE_LIMIT ){
                            objectFound = areabig;
                            *frame_processedROI=0;
                            *xposROI=xtemp;
                            *yposROI=ytemp;
                            emit send_position(*xposROI, *yposROI , objectFound);
                        }
                        else{
                            objectFound = 0;}
                        (*frame_processedROI)++;
                    }else{objectFound = areabig;
                        *frame_processedROI=1;
                        *xposROI=xtemp;
                        *yposROI=ytemp;
                        emit send_position(*xposROI, *yposROI , objectFound);}
                    // }
                }
                else {

                    (*frame_processedROI)++;
                    objectFound = 0;}
            }else{
                (*frame_processedROI)++;
                objectFound = 0;
            }
            if(condition){
                emit send_position_conditioning(*xposROI,*yposROI,objectFound);
            }
        }
        emit send_clock(true);             
}
void tracker::stop()
{
}
void tracker::conditioning(bool conditioning_activated){
    condition=conditioning_activated;
}
