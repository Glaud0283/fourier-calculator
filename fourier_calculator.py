import sympy as sp
from datetime import datetime
import sys
from tabulate import tabulate
import matplotlib.pyplot as plt
sys.setrecursionlimit(10000)

a1, a2, a3, a4, a5 = 2734, -6718, 5682, -1931, 228
t = sp.symbols('t', real=True)
F_1 = 100 * (a1 * t**5 + a2 * t**4 + a3 * t**3 + a4 * t**2 + a5 * t)
F_2 = -500
T = 1.5


def calculate_a0(F, T, a, b):
    a0 = (1 / T) * (sp.integrate(F, (t, 0, 1))) / 2
    return a0.subs(sp.pi, sp.N(sp.pi))

def calculate_an(F, T, a, b, n):
    an = (2 / T) * sp.integrate(F * sp.cos(2 * n * sp.pi * t / T), (t, a, b))
    return an.subs(sp.pi, sp.N(sp.pi))

def calculate_bn(F, T, a, b, n):
    bn = (2 / T) * sp.integrate(F * sp.sin(2 * n * sp.pi * t / T), (t, a, b))
    return bn.subs(sp.pi, sp.N(sp.pi))

def info(func):
    def wrapper(k):
        before = datetime.now()
        result = func(k)
        after = datetime.now()
        dur = after - before
        print("\n")
        print(f"fourier series is calculated for k = {k}. Time Passed: {dur.seconds} seconds")
        print("---")
        print("---")
        return result
    return wrapper
@info
def fourier(k):
    a0 = calculate_a0(F_1, T, 0, 1) + calculate_a0(F_2, T, 1, 1.5)
    fourier_series = a0 

    for n in range(1, k + 1):
        an_1 = calculate_an(F_1, T, 0, 1, n)
        an_2 = calculate_an(F_2, T, 1, 1.5, n)

        bn_1 = calculate_bn(F_1, T, 0, 1, n)
        bn_2 = calculate_bn(F_2, T, 1, 1.5, n)

        an = an_1 + an_2
        bn = bn_1 + bn_2
        fourier_series += an * sp.cos(2 * n * sp.pi * t / T) + bn * sp.sin(2 * n * sp.pi * t / T)
    return fourier_series.evalf()
    




def plot_fourier():
    f = fourier(1)
    plot = sp.plot(f, (t, 0, 3), title="Fourier Series", xlabel="Time [sec]", ylabel="F(t) [N]", ylim = (-1000, 1000), show = False)
    plot.save(f"{1}_iteration")
    print(f"F(t) = {f}")

    f = fourier(3)
    plot = sp.plot(f, (t, 0, 3), title="Fourier Series", xlabel="Time [sec]", ylabel="F(t) [N]", ylim = (-1000, 1000), show = False)
    plot.save(f"{3}_iteration")
    print(f"F(t) = {f}")

    f = fourier(10)
    plot = sp.plot(f, (t, 0, 3), title="Fourier Series", xlabel="Time [sec]", ylabel="F(t) [N]", ylim = (-1000, 1000), show = False)
    plot.save(f"{10}_iteration")
    print(f"F(t) = {f}")
    f = fourier(20)
    plot = sp.plot(f, (t, 0, 3), title="Fourier Series", xlabel="Time [sec]", ylabel="F(t) [N]", ylim = (-1000, 1000), show = False)
    plot.save(f"{100}_iteration")
    print(f"F(t) = {f}")
    f = fourier(100)
    plot = sp.plot(f, (t, 0, 3), title="Fourier Series", xlabel="Time [sec]", ylabel="F(t) [N]", ylim = (-1000, 1000), show = False)
    plot.save(f"{100}_iteration")
    print(f"F(t) = {f}")

def calculate_amplitude(k):
    data = []
    freqs = []
    amps = []
    col_headers = ["Frekans", "Genlik"]
    for n in range(1, k + 1):
        an_1 = calculate_an(F_1, T, 0, 1, n)
        an_2 = calculate_an(F_2, T, 1, 1.5, n)

        bn_1 = calculate_bn(F_1, T, 0, 1, n)
        bn_2 = calculate_bn(F_2, T, 1, 1.5, n)

        an = an_1 + an_2
        bn = bn_1 + bn_2
        amp = round(sp.sqrt(an**2 + bn**2), 3) 
        freq = round(n / T, 3)
        data. append((freq, amp))
        freqs.append(freq)
        amps.append(amp)
    print(tabulate(data, headers=col_headers, tablefmt="fancy_grid", showindex="always"))
    plt.figure(figsize=(8, 4))
    plt.stem(freqs, amps, linefmt='none', markerfmt='bo', basefmt="none ")
    plt.title('Frequency Spectrum')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()
    

def plot_displacement():
    x = 0.345 * sp.cos(13.2*t - 0.051)
    plot = sp.plot(x, (t, 0, 10), title="ϴ(t)", xlabel="Time [sec]", ylabel="ϴ(t) ", ylim = (-0.5, 0.5), show = False)
    plot.save(f"displacement_graph")

def plot_f_a():
    c = 270
    k = 3500
    x = 0.345 * sp.cos(13.2*t - 0.051)
    v = sp.diff(x)
    F_a = -1 * (fourier(20) - x * k + v * c)
    plot = sp.plot(F_a, (t, 0, 10), title="A Noktasındaki Tepki Kuvveti", xlabel="Time [sec]", ylabel="F(t) [N] ", ylim = (-2000, 2000), show = False)
    plot.save(f"reactionforce_graph")

plot_f_a()
