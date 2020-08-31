#!/usr/bin/env python

from __future__ import division

from os import system, name

from topy.utils import get_logger

logger = get_logger(__name__)

if name == 'posix':
    system('rm -f Q4.K')
    system('rm -f Q4bar.K')
    system('rm -f Q4T.K')
    system('rm -f Q5B.K')
    system('rm -f H8.K')
    system('rm -f H8T.K')
    system('rm -f H18B.K')
elif name == 'win32':
    system('del Q4.K')
    system('del Q4bar.K')
    system('del Q4T.K')
    system('del Q5B.K')
    system('del H8.K')
    system('del H8T.K')
    system('del H18B.K')

logger.info('This may take a while, perhaps a few minutes...')
system('python3 Q4_K.py')
logger.info('1 of 7 done!')
system('python3 Q4bar_K.py')
logger.info('2 of 7 done!')
system('python3 Q4T_K.py')
logger.info('3 of 7 done!')
system('python3 Q5B_K.py')
logger.info('4 of 7 done! All 2D matrices created. Now 3D...')
system('python3 H8_K.py')
logger.info('5 of 7 done!')
system('python3 H8T_K.py')
logger.info('6 of 7 done!')
system('python3 H18B_K.py')
logger.info('7 of 7 done! All 3D matrices created. Finished.')

# EOF recreate_all.py
