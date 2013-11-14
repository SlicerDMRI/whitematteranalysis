import pstats
print " = full dirs:"

p = pstats.Stats('profile.txt')
p.strip_dirs().sort_stats('time').print_stats(20)


p.strip_dirs().sort_stats('calls').print_stats(40)

