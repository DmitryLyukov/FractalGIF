#include <string>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <thread>
#include <atomic>

#include "CImg.h"


long double iter_for_point(const long double Re, const long double Im,
        const int MAX_ITER = 2500) {
    
    std::complex<long double> c(Re, Im);
    int i = 1;
    std::complex<long double> z_i = c;
    long double sqr_abs = real(z_i) * real(z_i) + imag(z_i) * imag(z_i);

    while ((i < MAX_ITER) && (sqr_abs < 4)) {
        z_i = z_i * z_i + c;
        ++i;
        sqr_abs = real(z_i) * real(z_i) + imag(z_i) * imag(z_i);
    }
    
    return static_cast<long double>(i) / MAX_ITER;
}

void find_colors(const long double norm_clr,
        unsigned char &red, unsigned char &green, unsigned char &blue) {
   
    red   = static_cast<unsigned char>(pow(sin(4.6 * norm_clr), 2) * 255);
    green = static_cast<unsigned char>(pow(sin(4.6 * norm_clr * norm_clr * norm_clr), 2) * 255);
    blue  = static_cast<unsigned char>(pow(sin(4.6 * norm_clr * norm_clr), 2) * 255);
}


void painting(const long double x1, const long double x2,
        const long double y1, const long double y2,
        const int32_t width, const int32_t height,
        const char* file_name) {
    
    cimg_library::CImg<unsigned char> img(width, height, 1, 3);
    
    const long double dy = (y2 - y1) / height;
    const long double dx = (x2 - x1) / width;
    
    long double y = y2;
    
    for (int32_t i = 0; i < height; ++i) {
        long double x = x1;
                
        for (int32_t j = 0; j < width; ++j) {
            const long double norm_clr = 1. - iter_for_point(x, y);
            
            unsigned char red, green, blue;
            find_colors(norm_clr, red, green, blue);
            
            img(j, i, 0) = red;
            img(j, i, 1) = green;
            img(j, i, 2) = blue;
            
            x += dx;
        }
        
        y -= dy;
    }
    
    img.save_png(file_name);
}

void painting_thread(std::atomic<int64_t> &atom_cnt) {
    const size_t PICT_CNT = 500;
    const long double ZOOM_CF = 0.95;
    while (atom_cnt < PICT_CNT) {
        int64_t j = atom_cnt.fetch_add(1);
        if (j >= PICT_CNT) {
            break;
        }
        const long double x = -0.55040561843894;
        const long double y = -0.62735213178491;
        
        const long double x1 = x - (x + 2.) * pow(ZOOM_CF, j);
        const long double x2 = x + (1. - x) * pow(ZOOM_CF, j);
        const long double y1 = y - (y + 0.84375) * pow(ZOOM_CF, j);
        const long double y2 = y + (0.84375 - y) * pow(ZOOM_CF, j);
        
        std::string str = std::to_string(j);
        while (str.length() < 5) {
            str = '0' + str;
        }
        str = str + '_' + std::to_string(x1) + '_' + std::to_string(x2)
                  + '_' + std::to_string(y1) + '_' + std::to_string(y2)
                  + ".png";
        painting(x1, x2, y1, y2, 1280, 720, str.c_str());
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
