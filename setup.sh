#!/bin/bash

function set_samweb() {
    exp=sbn
    export SAM_EXPERIMENT=$exp
    export SAM_GROUP=$exp
    export SAM_STATION=$exp
    export SAM_WEB_HOST=sam$exp.fnal.gov
    export IFDH_BASE_URI=http://sam$exp.fnal.gov:8480/sam/$exp/api/
    setup sam_web_client
}

function get_proxycert() {
    kx509
    voms-proxy-init -noregen -rfc -voms 'fermilab:/fermilab/sbnd/Role=Analysis'
}

function set_spine_sbnana() {
    mrb newDev -f -v v09_78_06 -q e20:prof
    source localProducts_larsoft_v09_78_06_e20_prof/setup
    cd srcs
    mrb g -b feature/jlarkin_2024_numu sbnana
    mrb g -b feature/jlarkin_2024_numu sbnanaobj
    cd sbnanaobj
    git cherry-pick ed3bf08865794763dd2cd92ba67dca919349e3d9
    cp /exp/sbnd/app/users/mueller/spineprod/cafana/srcs/sbnanaobj/sbnanaobj/StandardRecord/classes_def.xml sbnanaobj/StandardRecord/classes_def.xml
    cd ../..
    mrbsetenv
    mrb i -j4
}

function set_spine_anaplot() {
    git clone https://github.com/justinjmueller/spine_anaplot_tools.git
    cd spine_anaplot_tools
    git checkout feature/rlazur_sbnd2025
    cd ..
}


set_samweb
get_proxycert
#set_spine_sbnana
#set_spine_anaplot
source localProducts_larsoft_v09_78_06_e20_prof/setup
mrbsetenv
mrbslp
