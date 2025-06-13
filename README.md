# 🧨 UOC Option Pricing Engine
Welcome to the thrilling world of Up-and-Out Call Option Pricing, where we simulate financial instruments with the same elegance as a 10-year-old learning to rollerblade. This project uses finite difference methods under a jump-diffusion model to price exotic options that normal people don’t understand and even fewer should trade.

## 📦 Contents
```bash
src/
├── core/                    # Where your option logic lives, free of emotion
├── pricing/                # Tridiagonal solvers because diagonals aren't enough
├── utils/                  # Math that makes your high school teacher cry
main.cpp                    # Just enough to run something, barely
Makefile                   # You still have to type 'make', sorry
```

## 🛠 Requirements
- g++ (C++17, because we have standards now)

- GSL (libgsl-dev, install it like a grown-up)

- OpenMP (optional, for people who want results today)

## 🚀 Build It Like You Mean It
```bash
make
```

Try not to blink. It might just compile.

## 🧪 Try It Out
Run the binary. If you see a price, congratulations. If not, well... that's what debuggers are for.
```bash
./main
```
Or write your own tests. I believe in you, kind of.

## 💡 Features
- Choose between a home-brewed tridiagonal solver or the GSL version because variety is the spice of bugs.

- Handles jumps, diffusions, and your hopes and dreams.

- Mesh-based solution grid so dense it should probably file taxes.

- Helper math functions that require external libraries and maybe a prayer.

## 🧹 Clean Up
```bash
make clean
```

Deletes your hopes, dreams, and the bin/ directory.

## 📌 Notes
- This is for academic and recreational pain only.

- Expect performance to vary depending on whether your machine is powered by coffee or despair.

- Y ∈ [0, 1). Outside that, we throw runtime errors like candy.

## 🫠 Here's What I Regret
- Not wrapping GslTridiagonalSolver in a sanity-checking decorator. Trusting a black-box numerical library with your money is like handing your keys to a raccoon because it looked confident.

- Assuming Y will behave. Y is supposed to be in [0,1), but will silently betray you the moment you get cocky.

- Hardcoding grid sizes. I mean, yeah, it works, but adjusting the mesh should feel less like diffusing a bomb.

- Storing parameters in vectors like it’s 2004. params[3] is... wait, is that Y? Or is it theta? Don’t lie, you’ve guessed before.

- Calling init() in constructors like a nervous tic. One day, someone’s going to forget and everything’s going to default to $S = 100 because we didn’t believe in configuration files.

- Debug printouts as if we’re running a fax machine. dumpPrint() is cute until it starts gaslighting you.

- Not testing both solvers under duress. Sure, they work... until they don’t, and suddenly you’re pricing negative options.

- Using cmath without checking for domain errors. Mathematically undefined behavior is just code with extra vibes.
