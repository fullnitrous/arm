gcc -ferror-limit=100 -framework CoreVideo -framework IOKit -framework Cocoa -framework GLUT -framework OpenGL -I "src/raylib-5.0/src/" src/libraylib_macos.a src/main.c src/opw.c -lm -o bin/arm
