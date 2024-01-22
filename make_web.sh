emcc -s ASYNCIFY -s USE_GLFW=3 -I "src/raylib-5.0/src/" src/raylib-5.0/src/libraylib.a src/main.c src/opw.c -lm -o bin/arm.html
