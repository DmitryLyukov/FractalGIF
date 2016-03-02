#include <vector>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <string>

#include <thread>
#include <atomic>

#include "CImg.h"


double iter_for_point(const std::complex<double> &z, const std::complex<double> &c,
        const int MAX_ITER = 500) {
    
    int i = 0;
    std::complex<double> z_i = z;
    double sqr_abs = pow(real(z_i), 2) + pow(imag(z_i), 2);

    while ((i < MAX_ITER) && (sqr_abs < 4)) {
        z_i = z_i * z_i + c;
        ++i;
        sqr_abs = pow(real(z_i), 2) + pow(imag(z_i), 2);
    }
    
    return static_cast<double>(i) / MAX_ITER;
}


void find_colors(const double norm_clr,
        unsigned char &red, unsigned char &green, unsigned char &blue) {
   
    red   = static_cast<unsigned char>(pow(sin(4.6 * norm_clr), 2) * 255);
    green = static_cast<unsigned char>(pow(sin(4.6 * norm_clr * norm_clr * norm_clr), 2) * 255);
    blue  = static_cast<unsigned char>(pow(sin(4.6 * norm_clr * norm_clr), 2) * 255);
}


void painting(const double x1, const double x2,
        const double y1, const double y2,
        const int32_t width, const int32_t height,
        const char* file_name, const float re, const float im) {
    
    cimg_library::CImg<unsigned char> img(width, height, 1, 3);
    
    const double dy = (y2 - y1) / height;
    const double dx = (x2 - x1) / width;
       
    double y = y2;
    
    std::complex<double> c(re, im);
    
    for (int32_t i = 0; i < height / 2; ++i) {
        double x = x1;
        
        for (int32_t j = 0; j < width; ++j) {
            std::complex<double> z(x, y);
            const double norm_clr = 1. - iter_for_point(z, c);
            
            unsigned char red, green, blue;
            find_colors(norm_clr, red, green, blue);
            
            img(j, i, 0) = red;
            img(j, i, 1) = green;
            img(j, i, 2) = blue;
            
            img(width - j - 1, height - i - 1, 0) = red;
            img(width - j - 1, height - i - 1, 1) = green;
            img(width - j - 1, height - i - 1, 2) = blue;
            
            x += dx;
        }
        
        y -= dy;
    }
    
    img.save_png(file_name);
}


void painting_thread(std::atomic<int64_t> &atom_cnt) {
    const unsigned int num_of_pict = 200;
    while (atom_cnt < num_of_pict) {
        int64_t j = atom_cnt.fetch_add(1);
        if (j >= num_of_pict) {
            break;
        }
        
        const double dAlpha = 2 * 3.14159265358979323846 / num_of_pict;
        const double t      = j * dAlpha;
        const double re     = 0.5 * cos(t) - 0.25 * cos(2 * t);
        const double im     = 0.5 * sin(t) - 0.25 * sin(2 * t);
        
        std::string str = std::to_string(j);
        while (str.length() < 5) {
            str = '0' + str;
        }
        str = str + '_' + std::to_string(re) + '_' + std::to_string(im) + ".png";
        
        painting(-1.9, 1.9, -1.25, 1.25, 912, 600, str.c_str(), re, im);
        std::cout << j << std::endl;
    }
}


int main(int argc, char* argv[]) {
    
    std::atomic<int64_t> atom_cnt(0);
    
    std::thread thread1(painting_thread, std::ref(atom_cnt));
    std::thread thread2(painting_thread, std::ref(atom_cnt));
    std::thread thread3(painting_thread, std::ref(atom_cnt));
    std::thread thread4(painting_thread, std::ref(atom_cnt));

    thread1.join();
    thread2.join();
    thread3.join();
    thread4.join(); 
        
    return 0;
}
