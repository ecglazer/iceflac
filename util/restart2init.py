#!/usr/bin/env python
'''Usage: run this script on the same directory as flac output. It will modify
'time.rs' and '_contents.rs' so that the restart will be as from 0-th step.
'''

import numpy as np

print 'updating "_contents.save" to "_contents.rs"'
orig = open('_contents.save').readline().split()
orig[1] = '0'
orig[2] = '0.00'

open('_contents.rs', 'w').write('    '.join(orig) + '\n')


# Reset time
t = np.fromfile('time.rs', dtype=float)
sec2yr = 1.0 / (86400*365)
print 'original time was %g yr' % (t[0] * sec2yr)
print 'reset it to 0 yr'
t[0] = 0

t.tofile('time.rs')

print 'ready to restart flac'
