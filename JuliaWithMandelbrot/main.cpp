#include <vector>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
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


int iter_for_mand(const double x, const double y, const int MAX_ITER = 500) {
    
    std::complex<double> c(x,y);
    if ((x > -0.75) && (y < 0.65) && (x < 0.375) && (y > -0.65)) {
        const double abs_c = abs(c - 0.25);
        if (abs_c * 2 < 1 - (real(c) - 0.25) / abs_c) {
            return MAX_ITER;
        }
    }
    
    int i = 0;
    std::complex<double> z_i = 0;
    double sqr_abs = pow(real(z_i), 2) + pow(imag(z_i), 2);

    while ((i < MAX_ITER) && (sqr_abs < 4)) {
        std::complex<double> sqr = pow(z_i, 2);
        z_i = sqr + c;
        ++i;
        sqr_abs = pow(real(z_i), 2) + pow(imag(z_i), 2);
    }
    
    return i;
}


void find_colors(const double norm_clr,
        unsigned char &red, unsigned char &green, unsigned char &blue) {
   
    red   = static_cast<unsigned char>(pow(sin(4.6 * norm_clr), 2) * 255);
    green = static_cast<unsigned char>(pow(sin(4.6 * norm_clr * norm_clr * norm_clr), 2) * 255);
    blue  = static_cast<unsigned char>(pow(sin(4.6 * norm_clr * norm_clr), 2) * 255);
}


void painting(const double x1, const double x2,
        const double y1, const double y2,
        const size_t width, const size_t height,
        const char* file_name, const double re, const double im) {
    
    cimg_library::CImg<unsigned char> img(width, height, 1, 3);
    
    const double dy = (y2 - y1) / ((height - 2) / 2);
    const double dx = (x2 - x1) / width;
       
    double y = y2;
    
    std::complex<double> c(re, im);
    
    for (size_t i = 0; i < (height - 2) / 4; ++i) {
        double x = x1;
        
        for (size_t j = 0; j < width; ++j) {
            std::complex<double> z(x, y);
            const double norm_clr = 1. - iter_for_point(z, c);
            
            unsigned char red, green, blue;
            find_colors(norm_clr, red, green, blue);
            
            img(j, i, 0) = red;
            img(j, i, 1) = green;
            img(j, i, 2) = blue;
            
            img(width - j - 1, (height - 2) / 2 - i - 1, 0) = red;
            img(width - j - 1, (height - 2) / 2 - i - 1, 1) = green;
            img(width - j - 1, (height - 2) / 2 - i - 1, 2) = blue;
            
            x += dx;
        }
        
        y -= dy;
    }
    
    const unsigned char black[] = { 0,0,0 };
    img.draw_line(0, height / 2 - 1, width - 1, height / 2 - 1, black);
    img.draw_line(0, height / 2    , width - 1, height / 2    , black);
    
    y = 1.3;
    const double dy_m = 2.6 / ((height - 2) / 2);
    const double dx_m = 3.68 / width;
    
    for (size_t i = 0; i < (height - 2) / 4; ++i) {
        double x = -2.4;
        
        for (size_t j = 0; j < width; ++j) {
            const double norm_clr = 1 - iter_for_mand(x, y) / 500.;
            
            unsigned char red, green, blue;
            find_colors(norm_clr, red, green, blue);
            
            img(j, height / 2 + 1 + i, 0) = red;
            img(j, height / 2 + 1 + i, 1) = green;
            img(j, height / 2 + 1 + i, 2) = blue;
            
            img(j, height - i - 1, 0) = red;
            img(j, height - i - 1, 1) = green;
            img(j, height - i - 1, 2) = blue;
            
            x += dx_m;
        }
        
        y -= dy_m;
    }
    
    size_t a = static_cast<size_t>((re + 2.4) / 3.68 * width);
    size_t b = height / 2 + 1 + static_cast<size_t>((1.4 - im) / 2.8 * height / 2);
    
    const unsigned char red[] = { 255,0,0 };
    img.draw_circle(a, b, 5, red);
    
    img.save_png(file_name);
}


void painting_thread(std::atomic<int64_t> &atom_cnt) {
    const int num_of_pict = 200;
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
        
        painting(-2.13, 2.13, -1.5, 1.5, 1000, 1410, str.c_str(), re, im);
        std::cout << j << std::endl;
    }
}


int main(int argc, char* argv[]) {
    
    size_t num_of_threads;
    if (argc > 1) {
        num_of_threads = atol(argv[1]);
    } else {
        num_of_threads = std::thread::hardware_concurrency();
    }
    if (num_of_threads == 0) {
        num_of_threads = 1;
    }

    std::atomic<int64_t> atom_cnt(0);
    
    std::vector<std::thread> threads;
    for (size_t i = 0; i < num_of_threads; ++i) {
        threads.push_back(std::thread (painting_thread, std::ref(atom_cnt)));
    }
    
    for (auto& thr: threads) {
        thr.join();
    }
        
    return 0;
}
