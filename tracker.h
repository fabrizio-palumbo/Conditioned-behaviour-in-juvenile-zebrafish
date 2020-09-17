#ifndef TRACKER_H
#define TRACKER_H

#include <QtCore>
#include <QThread>
#include <opencv2/opencv.hpp>
#include "opencv2/highgui/highgui.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <opencv/highgui.h>
#include <opencv/cv.h>
#include <fstream>
#include <opencv/highgui.h>
#include <opencv/cv.h>
#include <stdio.h>
#include <dshow.h>
#include <Windows.h>// mandatory for I don't remember what...
#include "errors.h"

class tracker : public QObject{
    Q_OBJECT
public:
    explicit tracker(QObject *parent = 0): QObject(parent){}
    void init();
public slots:
    void start();
    void start_tracking( cv::Mat* resultROI,cv::Mat* threshROI,cv::Mat* frameROI, cv::Mat* backROI,cv::Mat* foreROI,cv::BackgroundSubtractor* bgROI,int* xposROI, int* yposROI,int* frame_processedROI);
    void stop();
    void conditioning (bool conditioning_activated);


signals:
    void send_position(int, int, int);
    void send_position_conditioning(int, int, int);
    void send_clock(bool);
};

#endif // TRACKER_H
