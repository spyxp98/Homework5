import matplotlib.pyplot as plt
import math
import time
from tqdm import *

func1 = lambda x: (x**2)+math.sqrt(x-1)
func1_int = lambda x: (x**3)/3 + (math.sqrt(x-1))**3 / (3/2)

func2 = lambda x: math.cos(21 * math.pi * x / 2)
func2_int = lambda x: 2 /(21 *  math.pi) * math.sin(21 * math.pi * x / 2)

def preload(filename):
   with open(f'{filename}', 'r') as f:
      temp = list()
      for line in f:
         temp_line = line.split()
         temp.append(float(temp_line[0]))
   return temp


def legendre_pol(n, x):
   if n == 0:
      return 1
   elif n == 1:
      return x
   else:
      return (((2*n - 1)/n) * x * legendre_pol(n-1, x) -
              ((n-1)/n) * legendre_pol(n-2, x))

def legendre_deriv(n, x):
   return (n/(x**2 - 1)) * (x * legendre_pol(n, x) - legendre_pol(n-1, x))

# k-th root of polinomial of order n
def legendre_root(k, n):
   return ((1-1/(8 * n**2) + 1/(8 * n **3)) 
            * math.cos(math.pi * (4 * k - 1)/(4 * n + 2)))



def legendre(func, func_int, eps):
   error = 10
   n = 0
   data = list()
   while error > eps:
      integral_sum = 0
      for i in range(n):
         x_i = legendre_root(i + 1, n)
         a_i = 2 / ((1 - (x_i)**2) * (legendre_deriv(n, x_i))**2)
         integral_sum += func(x_i) * a_i
      error = (abs((integral_sum - func_int(1) + func_int(-1))
                  /(func_int(1) - func_int(-1))))
      n += 5
      print(f'leg - {n} - {error}')
      data.append(error)
   return data
      
def chebyshev_pol(n, x):
   if n == 0:
      return 1
   elif n == 1:
      return x
   else:
      return 2 * x * chebyshev_pol(n-1, x) - chebyshev_pol(n-2, x)

def chebyshev_roots(n, k):
   return math.cos(math.pi * (k - 0.5) / n)

def chebyshev(func, func_int, eps):
   error = 10
   n = 1
   data = list()
   while error > eps:
      integral_sum = 0
      for i in range(n):
         x_i = chebyshev_roots(n, i + 1)
         # a_i = sum((1/(m + 1) - 1/(m - 1)) * chebyshev_pol(m, x_i) / (1 * n/2) for m in range(2, n - 1, 2))
         # a_i += 2 * chebyshev_pol(0, x_i) / n
         a_i = (math.pi / n) * math.sin(math.pi * (i+1 - 0.5)/n)
         integral_sum += a_i * func(x_i)
      error = (abs((integral_sum - func_int(1) + func_int(-1))
                  /(func_int(1) - func_int(-1))))
      n += 1
      print(f'cheb - {n} - {error}')
      data.append(error)
   return data



def trapez(func, func_int, a, b, eps):
   N = 1
   step = (b - a)/2
   data = list()
   integral_sum = 0
   error = 10
   while error > eps:
      integral_sum = 0
      curr_x = a
      N *= 2
      step = (b - a)/N
      while curr_x < b:
         try:
            integral_sum += step*(func(curr_x) + func(curr_x + step))/2
         except ZeroDivisionError:
            pass
         curr_x += step
      error = (abs((integral_sum - func_int(b) + func_int(a))
                  /(func_int(b) - func_int(a))))
      data.append(error)
   return data

def plot(filename):
   with open(f'{filename}', 'r') as f:
      data1 = list()
      data2 = list()
      data3 = list()
      epsilons = preload('eps')
      for line in f:
         temp = line.split(',')
         data1.append(float(temp[0]))
         data2.append(float(temp[1]))
         data3.append(float(temp[2]))
   plt.scatter(epsilons[:7], data1, color = 'red', marker = 'o')
   plt.xscale('log')
   plt.scatter(epsilons[:7], data2, color = 'green', marker = '^')
   plt.xscale('log')
   plt.scatter(epsilons[:7], data3, color = 'blue', marker = 'v')
   plt.xscale('log')
   plt.legend(['Trapezium', 'Legendre', 'Chebyshev'])
   plt.xlabel('relative error')
   plt.ylabel('time, miliseconds')
   plt.show()


def time_it(func, func_int):
   epsilons = preload('eps')
   with open('results', 'w') as f:
      # f.write('Трапеции, Лежандр, Чебышев \n')
      for i in tqdm(range(10)): 
         eps = epsilons[i]
         curr_time = time.time()
         data = trapez(func, func_int, -1, 1, eps)
         f.write(f'{(time.time() - curr_time)*1000},')
         curr_time = time.time()
         data = legendre(func, func_int, eps)
         f.write(f'{(time.time() - curr_time)*1000},')
         curr_time = time.time()
         data = chebyshev(func, func_int, eps)
         f.write(f'{(time.time() - curr_time)*1000}')
         f.write('\n')




if __name__ == '__main__':
   # data1 = trapez(func1, func1_int, 1, 5, eps)
   # data2 = legendre(func2, func2_int, eps)
   # print(data2)
   # plot(data1, data2)
   # time_it(func2, func2_int)
   plot('results')