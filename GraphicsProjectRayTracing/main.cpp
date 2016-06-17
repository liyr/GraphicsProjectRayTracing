#include <iostream>  
#include <opencv2/core/core.hpp>  
#include <opencv2/highgui/highgui.hpp>
#include "Space.h"
#include "kdtree.h"
#include "operation.h"
#include <random>
#include <ctime>
#include <thread>
#include "proces.h"

using namespace cv;




int main()
{
    //bg_mat = imread("./rainbow_texture679.jpg");
    marble_mat = imread("./marble_texture4665.jpg");
    wd_mat = imread("./wooden-floor-texture-wallpaper-10.jpg");
    sk_mat = imread("./free_space_galaxy_texture_by_lyshastra-d77gh54.jpg");
    eh_mat = imread("./earth.jpg");
    auto begin = clock();
	Mat img(length, width, CV_8UC3, Scalar::all(255));
	//randu(img, Scalar::all(0), Scalar::all(255));
    Space space(length, width, depth);

    preCalBrdf();


    while(1)
    {
        for (int i = 0; i < photonMapIteTimes; ++i)
        {
            LightTracing();
        }
        //TODO:save the map
        break;
        
    }

	//TODO: add objects(five sphere? + 1 bunny?) + Creature Texture
	//TODO: background light
	//TODO: (ray-tracing + photon mapping), multi-texture, soft shadows,Depth of field，mirroring, Anti-Aliasing 
    kdtree root_kdtree(&phoarray, 0, phoarray.end() - phoarray.begin(), 0);
    
    colorPaint(img, root_kdtree);
    std::cout << clock() - begin << std::endl;


	namedWindow("效果示例");
	imshow("效果示例", img); 
	waitKey(6000);
	imwrite(".\\ans13.png", img);

    system("pause");
	return 0;
}
