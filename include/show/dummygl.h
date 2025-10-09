#ifndef __SHOW_DUMMYGL_H__
#define __SHOW_DUMMYGL_H__

#define glColor3f(a, b, c) ((void)0)
#define glEnable(a) ((void)0)
#define glBindTexture(a, b) ((void)0)
#define glDisable(a) ((void)0)
#define glTexCoord1f(a) ((void)0)
#define glTexCoord2f(a, b) ((void)0)
#define glPixelStorei(a, b) ((void)0)
#define glTexParameteri(a, b, c) ((void)0)
#define glTexEnvf(a, b, c) ((void)0)
#define glTexEnvi(a, b, c) ((void)0)
#define glTexImage1D(a, b, c, d, e, f, g, h) ((void)0)
#define GLboolean int
#define glGetBooleanv(a, b) ((void)0)
#define GLubyte int
#define glColor3ubv(a) ((void)0)
#define glBegin(a) ((void)0)
#define glEnd() ((void)0)
#define glPointSize(a) ((void)0)
#define glLineWidth(a) ((void)0)
#define glVertex3f(a, b, c) ((void)0)
#define glVertex3d(a, b, c) ((void)0)
#define glColor4f(a, b, c, d) ((void)0)
#define glColor4d(a, b, c, d) ((void)0)
#define GLfloat int
#define GLint int
#define GLuint int
#define GLdouble double
#define GLenum int
#define glPushMatrix() ((void)0)
#define glPopMatrix() ((void)0)
#define glMultMatrix(a) ((void)0)
#define glMultMatrixd(a) ((void)0)
#define glLoadName(a) ((void)0)
#define GL_SELECT 0
#define glFlush() ((void)0)
#define glDrawBuffer(a) ((void)0)
#define glPolygonMode(a, b) ((void)0)
#define glLineStipple(a, b) ((void)0)
#define glMatrixMode(a) ((void)0)
#define glLoadIdentity() ((void)0)
#define glFinish() ((void)0)
#define glOrtho(a, b, c, d, e, f) ((void)0)
#define glBlendFunc(a, b) ((void)0)
#define glRasterPos3f(a, b, c) ((void)0)
#define glLogicOp(a) ((void)0)
#define glClearColor(a, b, c, d) ((void)0)
#define glClear(a) ((void)0)
#define glRotated(a, b, c, d) ((void)0)
#define glGetFloatv(a, b) ((void)0)
#define glTranslated(a, b, c) ((void)0)
#define GL_EXP 0
#define GL_EXP2 0
#define GL_LINEAR 0
#define glFogi(a, b) ((void)0)
#define glFogf(a, b) ((void)0)
#define glFogfv(a, b) ((void)0)
#define glHint(a, b) ((void)0)
#define glGetIntegerv(a, b) ((void)0)
#define GL_RENDER 0
#define glSelectBuffer(a, b) ((void)0)
#define glRenderMode(a) (0)
#define glInitNames() ((void)0)
#define glPushName(a) ((void)0)
#define glDepthFunc(a) ((void)0)
#define glViewport(a, b, c, d) ((void)0)
#define glLoadMatrixd(a) ((void)0)
#define glReadPixels(a, b, c, d, e, f, g) ((void)0)
#define glReadBuffer(a) ((void)0)
#define glFrustum(a, b, c, d, e, f) ((void)0)
#define glGetDoublev(a, b) ((void)0)
#define GL_BACK 0
#define GL_TEXTURE_ENV 0
#define GL_TEXTURE_ENV_MODE 0
#define GL_MODULATE 0
#define GL_UNSIGNED_BYTE 0
#define GL_RGBA 0
#define GL_TEXTURE_2D 0
#define glTexImage2D(a, b, c, d, e, f, g, h, i) ((void)0)
#define glGenTextures(a, b) ((void)0)

#define gluLookAt(a, b, c, d, e, f, g, h, i) ((void)0)
#define gluPerspective(a, b, c, d) ((void)0)
#define gluPickMatrix(a, b, c, d, e) ((void)0)
#define gluUnProject(a, b, c, d, e, f, g, h, i) ((void)0)
#define gluProject(a, b, c, d, e, f, g, h, i) ((void)0)
#define gluInit(a, b) ((void)0)

#define glutSwapBuffers() ((void)0)
#define glutGetWindow() (0)
#define glutPostRedisplay() ((void)0)
#define glutSetWindow(a) ((void)0)
#define glutGetModifiers() (0)
#define GLUT_ACTIVE_CTRL 0
#define GLUT_ACTIVE_ALT 0
#define GLUT_ACTIVE_SHIFT 0
#define glutCloseFunc(a) ((void)0)
#define glutMotionFunc(a) ((void)0)
#define glutKeyboardUpFunc(a) ((void)0)
#define glutDisplayFunc(a) ((void)0)
#define glutCreateWindow(a) (0)
#define glutInitWindowSize(a, b) ((void)0)
#define glutInitWindowPosition(a, b) ((void)0)
#define glutInitDisplayMode(a) ((void)0)
#define GLUT_DOWN 0
#define GLUT_LEFT_BUTTON 0
#define GLUT_MIDDLE_BUTTON 0
#define GLUT_RIGHT_BUTTON 0
#define glutFullScreen() ((void)0)
#define glutReshapeWindow(a, b) ((void)0)
#define glutInit(a, b) ((void)0)
#define glutPositionWindow(a, b) ((void)0)
#define glutTimerFunc(a, b, c) ((void)0)
#define glutMainLoop() ((void)0)
#define glutMainLoopEvent() ((void)0)

#define _glutBitmapString(a, b) ((void)0)

#define GLUI_WIDGET(a) class a { \
public: \
	void disable() {}; \
	void enable() {}; \
	const char *get_text() {return "";}; \
	void set_int_val(int) {}; \
	void set_speed(int) {}; \
	void set_alignment(int) {}; \
	void set_int_limits(int, int) {}; \
	void set_float_limits(float, float) {}; \
	void set_float_val(float) {}; \
	void reset() {}; \
};

#define GLUI_SPINNER_FLOAT 0
#define GLUI_EDITTEXT_TEXT 0
#define GLUI_ALIGN_LEFT 0
#define GLUI_ALIGN_RIGHT 0
#define GLUI_ALIGN_CENTER 0
#define GLUI_SPINNER_INT 0
#define GLUI_TRANSLATION_XY 0
#define GLUI_TRANSLATION_X 0
#define GLUI_TRANSLATION_Y 0
#define GLUI_TRANSLATION_Z 0

typedef void (*GLUI_Update_CB) (int id);

GLUI_WIDGET(GLUI_Button)
GLUI_WIDGET(GLUI_RadioButton)
GLUI_WIDGET(GLUI_RadioGroup)
GLUI_WIDGET(GLUI_Spinner)
GLUI_WIDGET(GLUI_Panel)
GLUI_WIDGET(GLUI_Listbox)
GLUI_WIDGET(GLUI_EditText)
GLUI_WIDGET(GLUI_Rotation)
GLUI_WIDGET(GLUI_Checkbox)

class GLUI {
public:
	int get_glut_window_id() { return 0; };
	GLUI_Panel* add_panel(const char *) { return nullptr; };
	void set_main_gfx_window(int) {};
	GLUI_Spinner* add_spinner_to_panel(GLUI_Panel*, const char*, int, void*) { return nullptr; };
	GLUI_Spinner* add_spinner_to_panel(GLUI_Panel*, const char*, int, void*, int, GLUI_Update_CB) { return nullptr; };
	void sync_live() {};
	GLUI_Button* add_button_to_panel(GLUI_Panel *, const char*, int, GLUI_Update_CB) { return nullptr; };
	void add_separator() {};
	GLUI_Checkbox* add_checkbox_to_panel(GLUI_Panel*, const char*, int*) { return nullptr; };
	GLUI_Checkbox* add_checkbox_to_panel(GLUI_Panel*, const char*, int*, int, GLUI_Update_CB) { return nullptr; };
	void add_separator_to_panel(GLUI_Panel*) {};
	void add_statictext_to_panel(GLUI_Panel*, const char*) {};
	GLUI_Button* add_button(const char*, int, GLUI_Update_CB) { return nullptr; };
	GLUI_Panel* add_rollout(const char*, bool) { return nullptr; };
	GLUI_EditText* add_edittext_to_panel(GLUI_Panel*, const char*, int, const char*) { return nullptr; };
	GLUI_Spinner* add_spinner(const char*, int, int*) { return nullptr; };
	GLUI_RadioButton* add_radiobutton_to_group(GLUI_RadioGroup*, const char*) { return nullptr; };
	GLUI_RadioGroup* add_radiogroup_to_panel(GLUI_Panel*, int*, int, GLUI_Update_CB) { return nullptr; };
	GLUI_RadioGroup* add_radiogroup_to_panel(GLUI_Panel*, int*) { return nullptr; };
	GLUI_Panel* add_rollout_to_panel(GLUI_Panel*, const char*, bool) { return nullptr; };
	void add_column_to_panel(GLUI_Panel*, bool) {};
	void add_translation_to_panel(GLUI_Panel*, const char*, int, float*, int, GLUI_Update_CB) {};
	GLUI_Rotation* add_rotation_to_panel(GLUI_Panel*, const char*, int*, int, GLUI_Update_CB) { return nullptr; };
	void add_checkbox(const char*, int*) {};
	void add_column(bool) {};
	void show() {};
};

class GLUI_Master {
public:
	void set_glutIdleFunc(void(&)()) {};
	void set_glutKeyboardFunc(void(&)(unsigned char, int, int)) {};
	void set_glutMouseFunc(void(&)(int, int, int, int)) {};
	void set_glutReshapeFunc(void(&)(int, int)) {};
	void set_glutSpecialFunc(void(&)(int, int, int)) {};
	GLUI* create_glui(const char *) { return nullptr; };
};
extern GLUI_Master GLUI_Master;

#endif
