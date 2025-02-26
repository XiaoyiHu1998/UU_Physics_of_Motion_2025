# Practical 2: Constraints

**Handout date:** 24-02-2025
**Deadline:** 11-03-2025 @ 08:59

## 1. Introduction

The second practical generalizes and extends the first practical by working with constraint-based velocity and position resolution. We still only work with rigid bodies. The objectives of the practical are:

1. Implement impulse-based velocity resolution by constraints (Topic 5).
    
2. Implement position correction by constraints (Topic 6).
    
3. Generalize collision resolution to work with the implemented constraint-resolution framework.
    
4. Extend the framework with some chosen effects.
    

Similarly to the first practical, this zip file contains the repository for the skeleton on which you will build your second practical. The environment is, for the most part, identical to the 1st practical. However, it is not advised to use your practical 1 implementation as a starting point for practical 2, since there are some important differences between the two templates.

## 2. Background

The practical is mostly about implementing procedures for velocity and position resolution to satisfy constraints. This resolution integrates into the game-physics engine as follows:

In each scene iteration:

1. Integrate the accelerations into velocities, and the velocities into positions and orientations (as in Practical 1)
    
2. Detect collisions
    
3. Resolve collisions through a generalized method that uses constraints, instead of the Practical 1 collision handling procedure.
    
4. Correct velocities to be tangent to user constraints (= constraints loaded from a file).
    
5. Correct positions to satisfy user constraints
    

The constraints we deal with in this practical are purely holonomic and bivariate. That is, each constraint is of the form , where and are two positions on two meshes (can be the same mesh). We distinguish between user constraints, read from a file, and collision constraints, created on-the-fly in each iteration. Moreover, we distinguish between equality constraints and inequality constraints . The difference between the two constraint types only matters for determining whether the constraint is violated, and otherwise should not affect velocities or positions. In practice, we use some tolerance that defines what we deem an 'acceptable' validity (opting for for equality constraints), rather than trying to achieve which is unattainable numerically.

### 2.1 Working with Rigid Bodies

Constraints can be defined between any two points on a body (which happen to be vertices in user constraints). Just as in practical 1, the bodies are rigid, and therefore the only degrees of freedom for a mesh are its COM position , orientation quaternion , linear COM velocity and angular velocity . Your implementation of constraint resolution should only work with and change these variables, and not touch any individual vertex. Specifically, never alter the vertex positions `currV` directly; these are updated for you by the template based on the current COM and orientation (as in practical 1).

### 2.2 Velocity Resolution

For equality constraints, the total velocities and should always satisfy , where is the gradient of the constraint, and is a vector comprising in order (note: this vector contains 12 variables). If , you will be computing to satisfy , using the Lagrange multiplier method (Topic 5). This requires setting up a (inverse) mass matrix, with the two objects’ masses and their (inverse) inertia tensors in order.

For fixed bodies, use for inverse mass and inverse inertia tensor. This will simulate the correct effect (as if they have infinite resistance). The part in the mass matrix corresponding to the linear velocity has the scalar masses and repeated 3 times each in the diagonal of the matrix, for the x, y, z components of the respective velocities.

### 2.3 Position Correction

Position correction is similar to velocity correction, except that we take the easy route (in the basic practical requirements), and only correct _linearly_. That is, we do not change , only of every body. This means that the mass matrix has a size of just and contains only body masses, without any inertia tensor components. The Jacobian only contains derivatives relating to linear movement.

## 3. Installation

This installation is exactly like that of Practical 1, repeated here for completeness.

### 3.1 Windows

On a Windows machine, use `cmake-gui` to compile the skeleton. Create a `build` folder inside the practical root directory. In `cmake-gui`, set the source code path and build directory, then configure and generate a Visual Studio solution. Open the solution and set the startup project to `Practical2_bin`.

### 3.2 MacOS / Linux

To compile the environment, run:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

## 4. The coding environment for the tasks

1. A constraint file is being read and has to be given as the third argument to the executable.
    
2. The function `handleCollision()` needs to be written by you to work with constraints.
    
3. `updateScene()` is already implemented.
    
4. Most work is in `Constraints.h`, where velocity and position correction functions need to be implemented.
    

## 5. Submission

Submit all modified files in a single zip file on Blackboard before **Tuesday, 11 March, 09:00 AM**. The practical will be checked during a session on the same day.

You are encouraged to work in pairs. Solo work requires prior permission.

## 6. Frequently Asked Questions

Q: I am getting "alignment" errors when compiling in Windows.  
A: Delete everything and reinstall using 64-bit configuration in `cmake-gui`.

Good luck! If you have questions, post them in the Practicals channel on Teams.