# Arm

Robotic arm engineering toolkit. For the theory behind how all this shit works
refer to the posts I have made on my website (still ongoing).

* [B-Splines and NURBS](https://fullnitrous.com/post/5RI2i)
* [Arcs and Bezier curves](https://fullnitrous.com/post/Ss9rw)
* [Spline Interpolation](https://fullnitrous.com/post/yyRW5)
	* [Visualizer](https://fullnitrous.com/post/gvj8k)
* [Forward Kinematics](https://fullnitrous.com/post/yKXjy)
	* [Visualizer](https://fullnitrous.com/post/EeEyd)
* [Inverse Kinematics](https://fullnitrous.com/post/qWlH6)

## Repository Structure

* `src/bin`
	* Binary files for compiling core source files and some libraries.
* `src/core`
	* Source code (ANSI C) for all the robotic arm specific stuff.
* `src/editor`
	* The robotic arm visualizer, called editor currently because that is
	going to be its final usage. As a editor for the robotic arm.
* `src/lib`
	* Statically linked libraries.
* `src/makefile`
	* Has targets: `all` for everything, `tests` `src/core`
	unit testing, `editor` for the editor and `editor_web` for the
	web version. Also `clean` for removing all generated files, includes
	files that aren't compiled by the makefile as well.
* `src/imgui-1.90.`, `raylib-5.0`, `rlimgui`
	* Dependencies (cloned repositories) for the editor.
* `src/tests`
	* Unit tests for `src/code`. Tests do not follow a specific C
	standard like `src/core` files do.
