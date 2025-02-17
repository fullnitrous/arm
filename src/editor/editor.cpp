#define RL_CULL_DISTANCE_NEAR 0.01
#define RL_CULL_DISTANCE_FAR  1000.0

#include <stdio.h>
#include <math.h>
#include <raylib.h>
#include <raymath.h>
#include <rcamera.h>
#include <rlgl.h>
#include <imgui.h>
#include <rlImGui.h>

#include <opw.h>
#include <util.h>
#include <intrpl.h>
#include <eval.h>

#define WINDOW_WIDTH  2880
#define WINDOW_HEIGHT 1620

#define COL_PURPLE (Color){127, 131, 255, 255}
#define COL_RED    (Color){255, 74,  152, 255}
#define COL_BLUE   (Color){10,  250, 250, 255}
#define COL_GRAY   (Color){127, 127, 127, 255}

#define COL_BG     (Color){24,   24,  24,  255}
#define COL_OBJ    (Color){48,   48,  48,  255}

#define CAMERA_MOVE_SPEED               0.09
#define CAMERA_ROTATION_SPEED           0.03
#define CAMERA_PAN_SPEED                0.2
#define CAMERA_MOUSE_MOVE_SENSITIVITY   0.003
#define CAMERA_MOUSE_SCROLL_SENSITIVITY 1.5
#define CAMERA_ORBITAL_SPEED            0.5

void update_camera(Camera *camera, int disable_keys, int disable_mouse) {
	Vector2 mouse_delta = GetMouseDelta();
	
	if(!disable_keys) {
		if(IsKeyDown(KEY_J)) { CameraPitch(camera, -CAMERA_ROTATION_SPEED, 0, 0, 0); }
		if(IsKeyDown(KEY_K)) { CameraPitch(camera, CAMERA_ROTATION_SPEED, 0, 0, 0); }
		if(IsKeyDown(KEY_L)) { CameraYaw(camera, -CAMERA_ROTATION_SPEED, 0); }
		if(IsKeyDown(KEY_H)) { CameraYaw(camera, CAMERA_ROTATION_SPEED, 0); }
		if(IsKeyDown(KEY_Q)) { CameraRoll(camera, -CAMERA_ROTATION_SPEED); }
		if(IsKeyDown(KEY_E)) { CameraRoll(camera, CAMERA_ROTATION_SPEED); }
		if(IsKeyDown(KEY_W)) { CameraMoveForward(camera, CAMERA_MOVE_SPEED, 0); }
		if(IsKeyDown(KEY_A)) { CameraMoveRight(camera, -CAMERA_MOVE_SPEED, 0); }
		if(IsKeyDown(KEY_S)) { CameraMoveForward(camera, -CAMERA_MOVE_SPEED, 0); }
		if(IsKeyDown(KEY_D)) { CameraMoveRight(camera, CAMERA_MOVE_SPEED, 0); }
		if(IsKeyDown(KEY_SPACE)) { CameraMoveUp(camera, CAMERA_MOVE_SPEED); }
		if(IsKeyDown(KEY_LEFT_CONTROL)) { CameraMoveUp(camera, -CAMERA_MOVE_SPEED); }
		if(IsKeyPressed(KEY_KP_SUBTRACT)) { CameraMoveToTarget(camera, 2.0f); }
		if(IsKeyPressed(KEY_KP_ADD))      { CameraMoveToTarget(camera, -2.0f); }
	}
	
	if(!disable_mouse) {
		if(IsKeyDown(KEY_LEFT_SHIFT) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
			if(mouse_delta.x > 0.0f) { CameraMoveRight(camera, CAMERA_PAN_SPEED, 0); }
			if(mouse_delta.x < 0.0f) { CameraMoveRight(camera, -CAMERA_PAN_SPEED, 0); }
			if(mouse_delta.y > 0.0f) { CameraMoveUp(camera, -CAMERA_PAN_SPEED); }
			if(mouse_delta.y < 0.0f) { CameraMoveUp(camera, CAMERA_PAN_SPEED); }
		} else if(IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
			CameraYaw(camera, -mouse_delta.x*CAMERA_MOUSE_MOVE_SENSITIVITY, 0);
			CameraPitch(camera, -mouse_delta.y*CAMERA_MOUSE_MOVE_SENSITIVITY, 0, 0, 0);
		}
		CameraMoveToTarget(camera, -GetMouseWheelMove());
	}
	
	return;
}

ImFont* font;

void setup_raylib(void) {
	SetConfigFlags(FLAG_WINDOW_RESIZABLE);
	SetConfigFlags(FLAG_MSAA_4X_HINT);
	//SetConfigFlags(FLAG_WINDOW_HIGHDPI);
	InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "arm");
	
	SetTargetFPS(144);

    rlImGuiBeginInitImGui();
    ImGui::StyleColorsDark();
	
	ImGuiIO& io = ImGui::GetIO();
	io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;


	font = io.Fonts->AddFontFromFileTTF("./editor/font.ttf", 34);

    rlImGuiEndInitImGui();
}

Vector3 from_vec3t(vec3_t v) {
	return (Vector3){(float)v.x, (float)v.y, (float)v.z};
}

void draw_axes(void) {
	float l = 2;
	Vector3 origin = {0, 0, 0};
	Vector3 x_axis = {l, 0, 0};
	Vector3 y_axis = {0, l, 0};
	Vector3 z_axis = {0, 0, l};
	DrawLine3D(origin, x_axis, COL_RED);
	DrawLine3D(origin, y_axis, COL_PURPLE);
	DrawLine3D(origin, z_axis, COL_BLUE);
	rlDrawRenderBatchActive();
	return;
}

void raylib_matrix(double* mat, float* raymat) {
	for(int i = 0; i < 16; i++) {
		raymat[i] = mat[4*(i%4)+i/4];
	}
	return;
}

double degrees(double radians) {
	return radians*180/M_PI;
}

void draw_arm_joint(double angle) {
	double r = 0.3;
	DrawCylinderWires((Vector3){0, 0, 0}, 0.2, 0.2, 0, 20, COL_PURPLE);
	DrawCylinder((Vector3){0, -0.3, 0}, 0.0, 0.1, 0.3, 20, COL_PURPLE);
	DrawLine3D((Vector3){0, -0.3, 0}, (Vector3){0.3, -0.3, 0}, COL_BLUE);
	for(int i = 0; i < 20; i++) {
		double a1 = -(i/20.0)*angle;
		double a2 = -((i+1)/20.0)*angle;
		vec3_t p1 = (vec3_t){r*cos(a1), -0.3, r*sin(a1)};
		vec3_t p2 = (vec3_t){r*cos(a2), -0.3, r*sin(a2)};
		DrawLine3D(from_vec3t(p1), from_vec3t(p2), COL_BLUE);
	}
	rlPushMatrix();
	rlRotatef(degrees(angle), 0.0, 1.0, 0.0);
	DrawLine3D((Vector3){0, -0.3, 0}, (Vector3){0.3, -0.3, 0}, COL_BLUE);
	rlPopMatrix();
	return;
}

void draw_arm(Model* linkages, kinstate_t* state, double* frames) {
	
	/* real 3d model */
	DrawModel(linkages[0], (Vector3){0, 0, 0}, 0.1f, COL_OBJ);
	DrawModelWires(linkages[0], (Vector3){0, 0, 0}, 0.1f, GRAY);
	
	rlPushMatrix();
	rlRotatef(degrees(state->linkage[0].theta)-90, 0.0, 0.0, 1.0);
	DrawModel(linkages[1], (Vector3){0, 0, 0}, 0.1f, COL_OBJ);
	DrawModelWires(linkages[1], (Vector3){0, 0, 0}, 0.1f, GRAY);
	
	rlPushMatrix();
	rlTranslatef(0, 0.81, 0.81);
	rlRotatef(270-degrees(state->linkage[1].theta) + 18, 1.0, 0.0, 0.0);
	rlTranslatef(0, -0.81, -0.81);
	DrawModel(linkages[2], (Vector3){0, 0, 0}, 0.1f, COL_OBJ);
	DrawModelWires(linkages[2], (Vector3){0, 0, 0}, 0.1f, GRAY);
	rlPopMatrix();
	rlPopMatrix();

	/* mathematical model */
	float t1[16], t2[16];
	rlPushMatrix();
	rlPushMatrix();
	rlRotatef(degrees(state->linkage[0].theta), 0.0, 0.0, 1.0);
	rlRotatef(90, 1.0, 0.0, 0.0);
	draw_arm_joint(-state->linkage[0].theta);
	rlPopMatrix();
	raylib_matrix(frames, t2);
	DrawLine3D((Vector3){0,0,0}, (Vector3){t2[12], t2[13], t2[14]}, COL_GRAY);

	DrawModel(linkages[0], (Vector3){0, 0, 0}, 0.1f, COL_OBJ);
	DrawModelWires(linkages[0], (Vector3){0, 0, 0}, 0.1f, GRAY);	

	for(int i = 0; i < 5; i++) {
		rlPushMatrix();
		raylib_matrix(frames+16*i, t1);
		rlMultMatrixf(t1);
		rlPushMatrix();
		rlRotatef(degrees(state->linkage[i+1].theta), 0.0, 0.0, 1.0);
		rlRotatef(90, 1.0, 0.0, 0.0);
		draw_arm_joint(-state->linkage[i+1].theta);
		rlPopMatrix();
		rlPopMatrix();
		raylib_matrix(frames+16*(i+1), t2);
		DrawLine3D((Vector3){t1[12], t1[13], t1[14]}, (Vector3){t2[12], t2[13], t2[14]}, COL_GRAY);
	}
	rlPopMatrix();	

	return;
}

vec3_t position_function1(ppoly_t ppoly, double t) {
	return eval_ppoly1(&ppoly, t).p;
}

vec3_t position_function3(ppoly_t ppoly, double t) {
	return eval_ppoly3(&ppoly, t).p;
}

vec3_t position_function5(ppoly_t ppoly, double t) {
	return eval_ppoly5(&ppoly, t).p;
}

vec3_t (*position_function)(ppoly_t, double) = &position_function3;

rotm_t orientation_function(double t) {
	return opw_euler(0, M_PI, 0);
}

void draw_path(ppoly_t ppoly) {
	for(int i = 1; i < 200; i++) {
		double t0 = range(0, 1, i-1, 200);
		double t1 = range(0, 1, i, 200);
		vec3_t p0 = position_function(ppoly, t0);
		vec3_t p1 = position_function(ppoly, t1);
		DrawLine3D(from_vec3t(p0), from_vec3t(p1), COL_GRAY);
	}
}

int ScaleToDPII(int value)
{
    return int(GetWindowScaleDPI().x * value);
}

class Window
{
public:
	bool open = false;

	RenderTexture ViewTexture;

	virtual void setup() = 0;
	virtual void shutdown() = 0;
	virtual void show() = 0;
	virtual void update() = 0;

	bool Focused = false;

	Rectangle ContentRect = { 0 };
};

class SceneView : public Window {
	public:
		Camera3D camera = { 0 };

		double sols[6*8];	
		double frames[6*4*4];
		int sol = 0;	
		double scale = 0.01;
		opw_t arm = {
			.c1 = 81*scale, .c2 = 550*scale, .c3 = 450*scale, .c4 = 100*scale,
			.a1 = 81*scale, .a2 = -70*scale,
			.b = 0*scale
		};
		
		kinstate_t state = opw_kinstate(arm);
		
		vec3_t target_pos;
		rotm_t target_rot;
		
		double t = 0;
		
		vec3_t points[5];
		
		vec3_t linear_mem0[EVAL_PPOLY1_MEMSIZE(5)];
		double linear_mem1[5];
		ppoly_t linear_;
		
		vec3_t cubic_mem0[EVAL_PPOLY3_MEMSIZE(5)];
		double cubic_mem1[5];
		ppoly_t cubic_;
		int pivot3_mem[EVAL_PPOLY3_PIVOT_MEMSIZE(5)];
		double matrix3_mem[EVAL_PPOLY3_EQSYS_MEMSIZE(5)];
		double b3_mem[EVAL_PPOLY3_B_MEMSIZE(5)];
		
		vec3_t quintic_mem0[EVAL_PPOLY5_MEMSIZE(5)];
		double quintic_mem1[5];
		ppoly_t quintic_;
		int pivot5_mem[EVAL_PPOLY5_PIVOT_MEMSIZE(5)];
		double matrix5_mem[EVAL_PPOLY5_EQSYS_MEMSIZE(5)];
		double b5_mem[EVAL_PPOLY5_B_MEMSIZE(5)];
		
		ppoly_t active;
		
		Model linkages[3];
		
		int sign = 1;

		void setup() override {
			ViewTexture = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());
			
			camera.position   = (Vector3){10, 10, 10};
			camera.target     = (Vector3){0, 0, 0};
			camera.up         = (Vector3){0, 0, 1};
			camera.projection = CAMERA_PERSPECTIVE;
			camera.fovy       = 45.0;
			
			Image img = GenImageChecked(ScaleToDPII(256), ScaleToDPII(256), ScaleToDPII(32), ScaleToDPII(32), DARKGRAY, WHITE);
			GridTexture = LoadTextureFromImage(img);
			UnloadImage(img);
			GenTextureMipmaps(&GridTexture);
			SetTextureWrap(GridTexture, TEXTURE_WRAP_CLAMP);
			
			points[0] = (vec3_t){4,   4,  1};
			points[1] = (vec3_t){-4,  4, -1};
			points[2] = (vec3_t){-4, -4,  1};
			points[3] = (vec3_t){4,  -4,  0};
			points[4] = (vec3_t){5,  -4,  1};
			
			linear_ = eval_ppoly_init(linear_mem0, linear_mem1, 5);
			cubic_ = eval_ppoly_init(cubic_mem0, cubic_mem1, 5);
			quintic_ = eval_ppoly_init(quintic_mem0, quintic_mem1, 5);
				
			luctx_t ctx3 = eval_ppoly_luctx_init(pivot3_mem, matrix3_mem, b3_mem);
			luctx_t ctx5 = eval_ppoly_luctx_init(pivot5_mem, matrix5_mem, b5_mem);
			
			intrpl_ppoly1(points, 5, &linear_);
			intrpl_ppoly3(points, 5, NULL, NULL, &ctx3, &cubic_); 
			intrpl_ppoly5(points, 5, NULL, NULL, NULL, NULL, &ctx5, &quintic_); 
			
			active = cubic_;
			
			linkages[0] = LoadModel("./editor/base.obj");
			linkages[1] = LoadModel("./editor/shoulder.obj");
			linkages[2] = LoadModel("./editor/upper.obj");	
		}

		void shutdown() override {
			UnloadRenderTexture(ViewTexture);
			UnloadTexture(GridTexture);
		}

		void show() override {
			ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
			ImGui::SetNextWindowSizeConstraints(ImVec2(ScaleToDPII(400), ScaleToDPII(400)), ImVec2((float)GetScreenWidth(), (float)GetScreenHeight()));

			if (ImGui::Begin("Robot", &open, ImGuiWindowFlags_NoScrollbar)) {
				Focused = ImGui::IsWindowFocused(ImGuiFocusedFlags_ChildWindows);
				rlImGuiImageRenderTextureFit(&ViewTexture, true);
			}

			ImGui::End();
			ImGui::PopStyleVar();
		}

		void update() override {
			t += sign*0.0005;
			if(t > 1.0) { t -= 0.002; sign *= -1; }
			if(t < 0.0) { t += 0.002; sign *= -1; }
				
			if(!open) return;
			
			if (IsWindowResized()) {
				UnloadRenderTexture(ViewTexture);
				ViewTexture = LoadRenderTexture(GetScreenWidth(), GetScreenHeight());
			}
			
			BeginTextureMode(ViewTexture);
			ClearBackground(COL_BG);
			
			update_camera(&camera, 0, 0);	
			
			BeginMode3D(camera);
			
			for(int i = 0; i < 5; i++) {
				DrawSphere(from_vec3t(points[i]), 0.1, COL_PURPLE);
			}
			
			draw_axes();
			
			KeyboardKey key = (KeyboardKey)GetKeyPressed();
			switch(key) {
				case KEY_ONE:   active = linear_; position_function = &position_function1; break;
				case KEY_TWO:   active = cubic_; position_function = &position_function3; break;
				case KEY_THREE: active = quintic_; position_function = &position_function5; break;
				default: break;
			}
			
			// main logic
			target_pos = position_function(active, t);
			target_rot = orientation_function(t);	
			draw_path(active);
			opw_inverse(arm, target_pos, target_rot, sols);	
			opw_forward(state, frames);
			opw_update_kinstate(&state, sols + 6*sol);
			draw_arm(linkages, &state, frames);						
			// main logic
			
			EndMode3D();
			EndTextureMode();
		}
		
		Texture2D GridTexture = { 0 };
};

int main(void) {
	setup_raylib();	
	
	SceneView sv;
	sv.setup();
	sv.open = true;

	while(!WindowShouldClose()) {	
		BeginDrawing();
		ClearBackground(COL_BG);
	
		static bool run = true;
		static bool showDemoWindow = true;
		
		rlImGuiBegin();
		ImGui::PushFont(font);
	
		ImGui::DockSpaceOverViewport(0,  NULL, ImGuiDockNodeFlags_PassthruCentralNode);

		if(ImGui::BeginMainMenuBar()) {
			if (ImGui::BeginMenu("File")) {
				if (ImGui::MenuItem("Quit"))
					run = false;
				ImGui::EndMenu();
            }

			if (ImGui::BeginMenu("Window")) {
                if (ImGui::MenuItem("Demo Window", nullptr, showDemoWindow))
					showDemoWindow = !showDemoWindow;
                ImGui::EndMenu();
            }
			ImGui::EndMainMenuBar();
		}

		if(showDemoWindow) {
			ImGui::ShowDemoWindow(&showDemoWindow);
		}

		if(ImGui::Begin("Motion Profile")) {
			ImGui::TextUnformatted("Another window");
			ImGui::End();
		}

		if(ImGui::Begin("Profiles")) {
			if (ImGui::BeginTabBar("MyTabBar", ImGuiTabBarFlags_None)) {
				if (ImGui::BeginTabItem("Projects")) {
					ImGui::Text("1");
					ImGui::EndTabItem();
				}
				if (ImGui::BeginTabItem("Broccoli")) {
					ImGui::Text("2");
					ImGui::EndTabItem();
				}
				if (ImGui::BeginTabItem("Cucumber")) {
					ImGui::Text("3");
					ImGui::EndTabItem();
				}
				ImGui::EndTabBar();
			}
			ImGui::TextUnformatted("Another window");
			ImGui::End();
		}


		sv.update();
		sv.show();

		ImGui::PopFont();
		rlImGuiEnd();
		
		EndDrawing();
	}

	CloseWindow();

	return 0;
}
