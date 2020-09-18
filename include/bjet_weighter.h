// Copyright [2020] Tyler Mitchell

#ifndef INCLUDE_BJET_WEIGHTER_H_
#define INCLUDE_BJET_WEIGHTER_H_

#include <math.h>

#include <string>

/*
    Only implement first row due to our selection
    ##################################################################
    Event weight matrix:            b-tagged jets
    ------------------------------------------------------------------
    nBJets                |    0        1             2
    ------------------------------------------------------------------
      0                   |    1      1-SF      (1-SF1)(1-SF2)
                          |
      1                   |    0       SF    SF1(1-SF2)+(1-SF1)SF2
                          |
      2                   |    0        0           SF1SF2
    ##################################################################
  */

enum bveto_wp { loose, medium };

class bjet_weighter {
   private:
    int year, working_point;
    double current_wp;
    double loose_wp_2016, loose_wp_2018;
    double medium_wp_2016, medium_wp_2018;
    double get_sf(double, double, std::string);
    double get_loose_sf_2016(double, double, std::string);
    double get_loose_sf_2018(double, double, std::string);
    double get_medium_sf_2016(double, double, std::string);
    double get_medium_sf_2018(double, double, std::string);

   public:
    bjet_weighter(int _year, int _wp);
    ~bjet_weighter() {}

    double find_weight(double, double, double, double, double, double, std::string);
};

bjet_weighter::bjet_weighter(int _year, int _wp)
    : year(_year), working_point(_wp), loose_wp_2016(0.2217), loose_wp_2018(0.1241), medium_wp_2016(0.6321), medium_wp_2018(0.4184) {
    if (working_point != bveto_wp::loose && working_point != bveto_wp::medium) {
        throw "working point must be bveto_wp::loose or bveto_wp::medium";
    }

    if (year == 2016) {
        if (working_point == bveto_wp::loose) {
            current_wp = loose_wp_2016;
        } else {
            current_wp = medium_wp_2016;
        }
    } else {
        if (working_point == bveto_wp::loose) {
            current_wp = loose_wp_2018;
        } else {
            current_wp = medium_wp_2018;
        }
    }
}

double bjet_weighter::find_weight(double pt1, double flv1, double bs1, double pt2, double flv2, double bs2, std::string syst = "none") {
    if (pt1 > 0 && bs1 > current_wp && pt2 > 0 && bs2 > current_wp) {
        return (1 - get_sf(pt1, flv1, syst)) * (1 - get_sf(pt2, flv2, syst));  // (1-SF1)*(1 - SF2)
    } else if (pt1 > 0 && bs1 > current_wp) {
        return 1 - get_sf(pt1, flv1, syst);  // 1 - SF
    } else if (pt2 > 0 && bs2 > current_wp) {
        return 1 - get_sf(pt2, flv2, syst);  // 1 - SF
    }
    return 1.;
}

double bjet_weighter::get_sf(double pt, double flv, std::string syst = "none") {
    if (working_point == bveto_wp::loose) {
        if (year == 2016) {
            return get_loose_sf_2016(pt, flv, syst);
        }
        // 2017 and 2018 are the same
        return get_loose_sf_2018(pt, flv, syst);
    } else if (working_point == bveto_wp::medium) {
        // medium working point
        if (year == 2016) {
            return get_medium_sf_2016(pt, flv, syst);
        }
        // 2017 and 2018 are the same
        return get_medium_sf_2018(pt, flv, syst);
    }
    return 1.;
}

double bjet_weighter::get_loose_sf_2016(double x, double flv, std::string syst = "none") {
    if (flv == 5) {
        auto nominal = 0.933791 * ((1. + (0.0115268 * x)) / (1. + (0.0103699 * x)));
        if (syst == "up" || syst == "down") {
            if (x < 30) {
                nominal += (syst == "up" ? +0.039079215377569199 : -0.039079215377569199);
            } else if (x < 50) {
                nominal += (syst == "up" ? +0.014356007799506187 : -0.014356007799506187);
            } else if (x < 70) {
                nominal += (syst == "up" ? +0.012236535549163818 : -0.012236535549163818);
            } else if (x < 100) {
                nominal += (syst == "up" ? +0.011882896535098553 : -0.011882896535098553);
            } else if (x < 140) {
                nominal += (syst == "up" ? +0.011785224080085754 : -0.011785224080085754);
            } else if (x < 200) {
                nominal += (syst == "up" ? +0.012215510942041874 : -0.012215510942041874);
            } else if (x < 300) {
                nominal += (syst == "up" ? +0.014544816687703133 : -0.014544816687703133);
            } else if (x < 600) {
                nominal += (syst == "up" ? +0.026683652773499489 : -0.026683652773499489);
            } else if (x < 1000) {
                nominal += (syst == "up" ? +0.055047694593667984 : -0.055047694593667984);
            }
        }
        return nominal;
    } else if (flv == 4) {
        auto nominal = 0.933791 * ((1. + (0.0115268 * x)) / (1. + (0.0103699 * x)));
        if (syst == "up" || syst == "down") {
            if (x < 30) {
                nominal += (syst == "up" ? +0.097698040306568146 : -0.097698040306568146);
            } else if (x < 50) {
                nominal += (syst == "up" ? +0.035890020430088043 : -0.035890020430088043);
            } else if (x < 70) {
                nominal += (syst == "up" ? +0.030591338872909546 : -0.030591338872909546);
            } else if (x < 100) {
                nominal += (syst == "up" ? +0.029707241803407669 : -0.029707241803407669);
            } else if (x < 140) {
                nominal += (syst == "up" ? +0.029463060200214386 : -0.029463060200214386);
            } else if (x < 200) {
                nominal += (syst == "up" ? +0.030538776889443398 : -0.030538776889443398);
            } else if (x < 300) {
                nominal += (syst == "up" ? +0.036362040787935257 : -0.036362040787935257);
            } else if (x < 600) {
                nominal += (syst == "up" ? +0.066709131002426147 : -0.066709131002426147);
            } else if (x < 1000) {
                nominal += (syst == "up" ? +0.13761924207210541 : -0.13761924207210541);
            }
        }
        return nominal;
    } else {
        auto nominal = 1.06337 + -0.000276004 * x + 1.25504e-06 * x * x + -8.9312e-10 * x * x * x;
        if (syst == "up" || syst == "down") {
            nominal *= (syst == "up" ? (1 + (0.0421943 + 5.30087e-05 * x + -6.87049e-08 * x * x))
                                     : (1 - (0.0421943 + 5.30087e-05 * x + -6.87049e-08 * x * x)));
        }
        return nominal;
    }
    return 1.;
}

double bjet_weighter::get_loose_sf_2018(double x, double flv, std::string syst = "none") {
    if (flv == 5) {
        auto nominal = 0.917829 + (0.00298278 * (log(x + 19) * (log(x + 18) * (3 - (0.422392 * log(x + 18))))));
        if (syst == "up" || syst == "down") {
            if (x < 30) {
                nominal += (syst == "up" ? +0.062023099511861801 : -0.062023099511861801);
            } else if (x < 50) {
                nominal += (syst == "up" ? +0.013962121680378914 : -0.013962121680378914);
            } else if (x < 70) {
                nominal += (syst == "up" ? +0.013880428858101368 : -0.013880428858101368);
            } else if (x < 100) {
                nominal += (syst == "up" ? +0.013638468459248543 : -0.013638468459248543);
            } else if (x < 140) {
                nominal += (syst == "up" ? +0.011050660163164139 : -0.011050660163164139);
            } else if (x < 200) {
                nominal += (syst == "up" ? +0.011366868391633034 : -0.011366868391633034);
            } else if (x < 300) {
                nominal += (syst == "up" ? +0.011010468937456608 : -0.011010468937456608);
            } else if (x < 600) {
                nominal += (syst == "up" ? +0.037737511098384857 : -0.037737511098384857);
            } else if (x < 1000) {
                nominal += (syst == "up" ? +0.069150865077972412 : -0.069150865077972412);
            }
        }
        return nominal;
    } else if (flv == 4) {
        auto nominal = 0.917829 + (0.00298278 * (log(x + 19) * (log(x + 18) * (3 - (0.422392 * log(x + 18))))));
        if (syst == "up" || syst == "down") {
            if (x < 30) {
                nominal += (syst == "up" ? +0.15505774319171906 : -0.15505774319171906);
            } else if (x < 50) {
                nominal += (syst == "up" ? +0.03490530326962471 : -0.03490530326962471);
            } else if (x < 70) {
                nominal += (syst == "up" ? +0.034701071679592133 : -0.034701071679592133);
            } else if (x < 100) {
                nominal += (syst == "up" ? +0.034096170216798782 : -0.034096170216798782);
            } else if (x < 140) {
                nominal += (syst == "up" ? +0.027626650407910347 : -0.027626650407910347);
            } else if (x < 200) {
                nominal += (syst == "up" ? +0.02841717004776001 : -0.02841717004776001);
            } else if (x < 300) {
                nominal += (syst == "up" ? +0.027526171877980232 : -0.027526171877980232);
            } else if (x < 600) {
                nominal += (syst == "up" ? +0.094343781471252441 : -0.094343781471252441);
            } else if (x < 1000) {
                nominal += (syst == "up" ? +0.17287716269493103 : -0.17287716269493103);
            }
        }
        return nominal;
    } else {
        auto nominal = 1.41852 + -0.00040383 * x + 2.89389e-07 * x * x + -3.55101e-11 * x * x * x;
        if (syst == "up" || syst == "down") {
            nominal *= (syst == "up" ? (1 + (0.0559259 + 1.96455e-05 * x + -3.60571e-08 * x * x))
                                     : (1 - (0.0559259 + 1.96455e-05 * x + -3.60571e-08 * x * x)));
        }
        return nominal;
    }
    return 1.;
}

double bjet_weighter::get_medium_sf_2016(double x, double flv, std::string syst = "none") {
    if (flv == 5) {
        auto nominal = 0.653526 * ((1. + (0.220245 * x)) / (1. + (0.14383 * x)));
        if (syst == "up" || syst == "down") {
            if (x < 30) {
                nominal += (syst == "up" ? +0.043795019388198853 : -0.043795019388198853);
            } else if (x < 50) {
                nominal += (syst == "up" ? +0.015845479443669319 : -0.015845479443669319);
            } else if (x < 70) {
                nominal += (syst == "up" ? +0.014174085110425949 : -0.014174085110425949);
            } else if (x < 100) {
                nominal += (syst == "up" ? +0.013200919143855572 : -0.013200919143855572);
            } else if (x < 140) {
                nominal += (syst == "up" ? +0.012912030331790447 : -0.012912030331790447);
            } else if (x < 200) {
                nominal += (syst == "up" ? +0.019475525245070457 : -0.019475525245070457);
            } else if (x < 300) {
                nominal += (syst == "up" ? +0.01628459244966507 : -0.01628459244966507);
            } else if (x < 600) {
                nominal += (syst == "up" ? +0.034840557724237442 : -0.034840557724237442);
            } else if (x < 1000) {
                nominal += (syst == "up" ? +0.049875054508447647 : -0.049875054508447647);
            }
        }
        return nominal;
    } else if (flv == 4) {
        auto nominal = 0.653526 * ((1. + (0.220245 * x)) / (1. + (0.14383 * x)));
        if (syst == "up" || syst == "down") {
            if (x < 30) {
                nominal += (syst == "up" ? +0.13138505816459656 : -0.13138505816459656);
            } else if (x < 50) {
                nominal += (syst == "up" ? +0.047536440193653107 : -0.047536440193653107);
            } else if (x < 70) {
                nominal += (syst == "up" ? +0.042522255331277847 : -0.042522255331277847);
            } else if (x < 100) {
                nominal += (syst == "up" ? +0.039602756500244141 : -0.039602756500244141);
            } else if (x < 140) {
                nominal += (syst == "up" ? +0.038736090064048767 : -0.038736090064048767);
            } else if (x < 200) {
                nominal += (syst == "up" ? +0.058426573872566223 : -0.058426573872566223);
            } else if (x < 300) {
                nominal += (syst == "up" ? +0.048853777348995209 : -0.048853777348995209);
            } else if (x < 600) {
                nominal += (syst == "up" ? +0.10452167689800262 : -0.10452167689800262);
            } else if (x < 1000) {
                nominal += (syst == "up" ? +0.14962516725063324 : -0.14962516725063324);
            }
        }
        return nominal;
    } else {
        auto nominal = 1.09286 + -0.00052597 * x + 1.88225e-06 * x * x + -1.27417e-09 * x * x * x;
        if (syst == "up" || syst == "down") {
            nominal *= (syst == "up" ? (1 + (0.101915 + 0.000192134 * x + -1.94974e-07 * x * x))
                                     : (1 - (0.101915 + 0.000192134 * x + -1.94974e-07 * x * x)));
        }
        return nominal;
    }
    return 1.;
}

double bjet_weighter::get_medium_sf_2018(double x, double flv, std::string syst = "none") {
    if (flv == 5) {
        auto nominal = 0.909339 + (0.00354 * (log(x + 19) * (log(x + 18) * (3 - (0.471623 * log(x + 18))))));
        if (syst == "up" || syst == "down") {
            if (x < 30) {
                nominal += (syst == "up" ? +0.065904870629310608 : -0.065904870629310608);
            } else if (x < 50) {
                nominal += (syst == "up" ? +0.015055687166750431 : -0.015055687166750431);
            } else if (x < 70) {
                nominal += (syst == "up" ? +0.013506759889423847 : -0.013506759889423847);
            } else if (x < 100) {
                nominal += (syst == "up" ? +0.015106724575161934 : -0.015106724575161934);
            } else if (x < 140) {
                nominal += (syst == "up" ? +0.014620178379118443 : -0.014620178379118443);
            } else if (x < 200) {
                nominal += (syst == "up" ? +0.012161554768681526 : -0.012161554768681526);
            } else if (x < 300) {
                nominal += (syst == "up" ? +0.016239689663052559 : -0.016239689663052559);
            } else if (x < 600) {
                nominal += (syst == "up" ? +0.039990410208702087 : -0.039990410208702087);
            } else if (x < 1000) {
                nominal += (syst == "up" ? +0.068454340100288391 : -0.068454340100288391);
            }
        }
        return nominal;
    } else if (flv == 4) {
        auto nominal = 0.909339 + (0.00354 * (log(x + 19) * (log(x + 18) * (3 - (0.471623 * log(x + 18))))));
        if (syst == "up" || syst == "down") {
            if (x < 30) {
                nominal += (syst == "up" ? +0.19771461188793182 : -0.19771461188793182);
            } else if (x < 50) {
                nominal += (syst == "up" ? +0.045167062431573868 : -0.045167062431573868);
            } else if (x < 70) {
                nominal += (syst == "up" ? +0.040520280599594116 : -0.040520280599594116);
            } else if (x < 100) {
                nominal += (syst == "up" ? +0.045320175588130951 : -0.045320175588130951);
            } else if (x < 140) {
                nominal += (syst == "up" ? +0.043860536068677902 : -0.043860536068677902);
            } else if (x < 200) {
                nominal += (syst == "up" ? +0.036484666168689728 : -0.036484666168689728);
            } else if (x < 300) {
                nominal += (syst == "up" ? +0.048719070851802826 : -0.048719070851802826);
            } else if (x < 600) {
                nominal += (syst == "up" ? +0.11997123062610626 : -0.11997123062610626);
            } else if (x < 1000) {
                nominal += (syst == "up" ? +0.20536302030086517 : -0.20536302030086517);
            }
        }
        return nominal;
    } else {
        auto nominal = 1.6329 + -0.00160255 * x + 1.9899e-06 * x * x + -6.72613e-10 * x * x * x;
        if (syst == "up" || syst == "down") {
            nominal *= (syst == "up" ? (1 + (0.122811 + 0.000162564 * x + -1.66422e-07 * x * x))
                                     : (1 - (0.122811 + 0.000162564 * x + -1.66422e-07 * x * x)));
        }
        return nominal;
    }
    return 1.;
}

#endif  // INCLUDE_BJET_WEIGHTER_H_
