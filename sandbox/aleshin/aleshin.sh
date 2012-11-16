#!/bin/sh
(echo $*; cat) | ssh lbartho@obssf3.unige.ch env LD_LIBRARY_PATH=/unige/nag/SunOS_5.8_sun4u/lib:/unige/workshop_6.0.2/SUNWspro/lib aleshin/aleshin
