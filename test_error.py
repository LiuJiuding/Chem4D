a=[1,2,3]
try:
    print(a[3])
except IndexError:
    print(a[1])