#include <iostream>
#include "Mesh.h"
#include "GLProjector.h"
#include "Deformer.h"
#include "Parameterization.h"
#include <commdlg.h>
#include <string>
#include <opencv2/opencv.hpp>
// the way the mesh is rendered
enum EnumDisplayMode {
    WIREFRAME,
    HIDDENLINE,
    FLATSHADED,
    SMOOTHSHADED,
    COLORSMOOTHSHADED, 
	PARAMETERSHADED, 
    TEXTRURESHADED,
    SKELETONS
};

// variables
int displayMode = SKELETONS;                     // current display mode
int mainMenu, displayMenu,filename;                 // glut menu handlers
int winWidth, winHeight;                            // window width and height
double winAspect;                                   // winWidth / winHeight;
int lastX, lastY;                                   // last mouse motion position
bool leftDown, middleDown, shiftDown;               // mouse down and shift down flags
double sphi = 90.0, stheta = 45.0, sdepth = 10;     // for simple trackball
double xpan = 0.0, ypan = 0.0;                      // for simple trackball
double zNear = 1.0, zFar = 100.0;                   // clipping
double g_fov = 45.0;
Eigen::Vector3d g_center;
double g_sdepth;
Mesh mesh; // our mesh
std::vector<Eigen::Vector3d> origianlPosition;
VertexList originalvertex;
GLuint m_textureID=0;
cv::Mat m_texture;
#include <opencv2/highgui/highgui.hpp>  
// editing mode
enum Mode {
    Viewing,
    Selection,
    Moving
};

Mode currentMode = Viewing;
int downX, downY; // mouse down position
int selectedHandleIndex = -1; // the index of the handle
Deformer* deformer = NULL; // Deformer, make a natural deformation
Parameterization * m_parameterization =nullptr;
bool showparameterizationresult = false;
// functions
void SetBoundaryBox(const Eigen::Vector3d& bmin, const Eigen::Vector3d& bmax);
void InitGL();
void InitMenu();
void InitGeometry();
std::string GetFileName()
{
	OPENFILENAME ofn;       
	TCHAR szFile[MAX_PATH];         
							  
	ZeroMemory(&ofn, sizeof(OPENFILENAME));
	ofn.lStructSize = sizeof(OPENFILENAME);
	ofn.hwndOwner = NULL;
	ofn.lpstrFile = szFile;
	ofn.lpstrFile[0] = '\0';
	ofn.nMaxFile = sizeof(szFile);
	ofn.lpstrFilter = "All(*.*)\0*.*\0Text(*.txt)\0*.TXT\0\0";
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
	//ofn.lpTemplateName =  MAKEINTRESOURCE(ID_TEMP_DIALOG);  

	if (GetOpenFileName(&ofn))
	{
		OutputDebugString(szFile);
		OutputDebugString("\r\n");
	}
	return szFile;
}

// window related 
void MenuCallback(int value);
void ReshapeFunc(int width, int height);

// rendering functions
void DisplayFunc();
void DrawWireframe();
void DrawHiddenLine();
void DrawFlatShaded();
void DrawSmoothShaded();
void DrawColorSmoothShaded();
void DrawSelectionRect();
void DrawSelectedVertices();
void DrawParameterResult();
void DrawTexturesShaded();
void DrawBoundaryLoops();
void DrawSkeletons();
// input related glut functions
void KeyboardFunc(unsigned char ch, int x, int y);
void MouseFunc(int button, int state, int x, int y);
void MotionFunc(int x, int y);

void SelectVertexByRect(); // select ROI for editing
int StartMove(); // a single editing operation


void SetBoundaryBox(const Eigen::Vector3d& bmin, const Eigen::Vector3d& bmax) {
    double PI = 3.14159265358979323846;
    double radius = (bmax - bmin).norm();
    g_center = 0.5 * (bmin + bmax);
    zNear = 0.2 * radius / sin(0.5 * g_fov * PI / 180.0);
    zFar = zNear + 2.0 * radius;
    g_sdepth = zNear + radius;
    zNear *= 0.1;
    zFar *= 10;
    sdepth = g_sdepth;
}

// init openGL environment
void InitGL() {
    GLfloat light0Position[] = {0, 1, 0, 1.0};

    // initialize GLUT stuffs
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(500, 500);
    glutCreateWindow("Comp5411 Mesh Viewer");

    //glClearColor(0.0, 0.0, 0.0, 0.0);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glPolygonOffset(1.0, 1.0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_DIFFUSE);
    glLightfv(GL_LIGHT0, GL_POSITION, light0Position);
    glEnable(GL_LIGHT0);

    glutReshapeFunc(ReshapeFunc);
    glutDisplayFunc(DisplayFunc);
    glutKeyboardFunc(KeyboardFunc);
    glutMouseFunc(MouseFunc);
    glutMotionFunc(MotionFunc);
}

// init right-click menu
void InitMenu() {
    displayMenu = glutCreateMenu(MenuCallback);
    glutAddMenuEntry("Wireframe", WIREFRAME);
    glutAddMenuEntry("Hidden Line", HIDDENLINE);
    glutAddMenuEntry("Flat Shaded", FLATSHADED);
    glutAddMenuEntry("Smooth Shaded", SMOOTHSHADED);
    glutAddMenuEntry("Color Smooth Shaded", COLORSMOOTHSHADED);
	glutAddMenuEntry("Parameter Shaded", PARAMETERSHADED);
	glutAddMenuEntry("Texture Shaded", TEXTRURESHADED);

    mainMenu = glutCreateMenu(MenuCallback);
    glutAddSubMenu("Display", displayMenu);
    glutAddMenuEntry("Build Deformer", 101);
    glutAddMenuEntry("Exit", 99);
	
    glutAttachMenu(GLUT_RIGHT_BUTTON);

}

// init geometry (if no input argument is provided)
void InitGeometry() {
    const int VSIZE = 4;
    const int HESIZE = 12;
    const int FSIZE = 4;
    int i;
    Vertex* v[VSIZE];
    HEdge* he[HESIZE];
    Face* f[FSIZE];

    for (i = 0; i < VSIZE; i++) {
        v[i] = new Vertex();
        mesh.vList.push_back(v[i]);
    }
    v[0]->SetPosition(Eigen::Vector3d(0.0, 0.0, 0.0));
    v[1]->SetPosition(Eigen::Vector3d(10.0, 0.0, 0.0));
    v[2]->SetPosition(Eigen::Vector3d(0.0, 10.0, 0.0));
    v[3]->SetPosition(Eigen::Vector3d(0.0, 0.0, 10.0));

    v[0]->SetNormal(Eigen::Vector3d(-0.577, -0.577, -0.577));
    v[1]->SetNormal(Eigen::Vector3d(0.0, -0.7, -0.7));
    v[2]->SetNormal(Eigen::Vector3d(-0.7, 0.0, -0.7));
    v[3]->SetNormal(Eigen::Vector3d(-0.7, -0.7, 0.0));

    for (i = 0; i < FSIZE; i++) {
        f[i] = new Face();
        mesh.fList.push_back(f[i]);
    }

    for (i = 0; i < HESIZE; i++) {
        he[i] = new HEdge();
        mesh.heList.push_back(he[i]);
    }
    for (i = 0; i < FSIZE; i++) {
        int base = i * 3;
        Mesh::SetPrevNext(he[base], he[base + 1]);
        Mesh::SetPrevNext(he[base + 1], he[base + 2]);
        Mesh::SetPrevNext(he[base + 2], he[base]);
        Mesh::SetFace(f[i], he[base]);
    }
    Mesh::SetTwin(he[0], he[4]);
    Mesh::SetTwin(he[1], he[7]);
    Mesh::SetTwin(he[2], he[10]);
    Mesh::SetTwin(he[3], he[8]);
    Mesh::SetTwin(he[5], he[9]);
    Mesh::SetTwin(he[6], he[11]);
    he[0]->SetStart(v[1]);
    he[1]->SetStart(v[2]);
    he[2]->SetStart(v[3]);
    he[3]->SetStart(v[0]);
    he[4]->SetStart(v[2]);
    he[5]->SetStart(v[1]);
    he[6]->SetStart(v[0]);
    he[7]->SetStart(v[3]);
    he[8]->SetStart(v[2]);
    he[9]->SetStart(v[0]);
    he[10]->SetStart(v[1]);
    he[11]->SetStart(v[3]);
    v[0]->SetHalfEdge(he[3]);
    v[1]->SetHalfEdge(he[0]);
    v[2]->SetHalfEdge(he[1]);
    v[3]->SetHalfEdge(he[2]);
}


// GLUT menu callback function
void MenuCallback(int value) {
    switch (value) {
    case 99:
        exit(0);
        break;
    case 101:
        deformer = new Deformer(&mesh);
        std::cout << "Deformer Building Finished\n";
        break;
    default:
        displayMode = value;
        glutPostRedisplay();
        break;
    }
}


// GLUT reshape callback function
void ReshapeFunc(int width, int height) {
    winWidth = width;
    winHeight = height;
    winAspect = (double)width / (double)height;
    glViewport(0, 0, width, height);
    glutPostRedisplay();
}


// GLUT display callback function
void DisplayFunc() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(g_fov, winAspect, zNear, zFar);
    //glOrtho(-2.0, 2.0, -2.0, 2.0, zNear, zFar);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(xpan, ypan, -sdepth);
    glRotatef(-stheta, 1.0, 0.0, 0.0);
    glRotatef(sphi, 0.0, 1.0, 0.0);
    glTranslatef(-g_center[0], -g_center[1], -g_center[2]);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    switch (displayMode) {
    case WIREFRAME:
        DrawWireframe();        
        break;
    case HIDDENLINE:
        DrawHiddenLine();
        break;
    case FLATSHADED:
        DrawFlatShaded();
        break;
    case SMOOTHSHADED:
        DrawSmoothShaded();
        break;
    case COLORSMOOTHSHADED:
        DrawColorSmoothShaded();
        break;
	case PARAMETERSHADED:
		DrawParameterResult();
		break;	
	case TEXTRURESHADED:
		DrawTexturesShaded();
		break;
    case SKELETONS:
        DrawSkeletons();
        break;
    }
    DrawBoundaryLoops();
    DrawSelectedVertices();
    if (currentMode == Selection && downX != lastX && downX != lastY) {
        DrawSelectionRect();
    }

    glutSwapBuffers();
}

//boundaryloop render function
//use different color to draw boundaryloops
void DrawBoundaryLoops() {
    HEdgeLGroup BoundaryLoop = mesh.BoundaryLoop;
    HEdgeList* abc;
    glBegin(GL_LINES);
    size_t i;
    size_t k;
    //glColor3f(1.0f, 0.0f, 1.0f);
    //for (std::vector<HEdgeList*>::const_iterator iter = BoundaryLoop.begin(); iter != BoundaryLoop.end();iter++) {
    for (k = 0; k < BoundaryLoop.size(); k++) {
        abc = BoundaryLoop[k];
        switch (k) {
            case(0):
                glColor3f(1.0f, 0.0f, 0.0f); break;
            case(1):
                glColor3f(0.0f, 1.0f, 0.0f); break;
            case(2):
                glColor3f(0.0f, 0.0f, 1.0f); break;
            case(3):
                glColor3f(1.0f, 1.0f, 0.0f); break;
            case(4):
                glColor3f(1.0f, 0.0f, 1.0f); break;
            case(5):
                glColor3f(0.0f, 1.0f, 1.0f); break;
            case(6):
                glColor3f(1.0f, 0.5f, 0.0f); break;
            case(7):
                glColor3f(1.0f, 0.0f, 0.5f); break;
            case(8):
                glColor3f(0.5f, 1.0f, 0.0f); break;
            case(9):
                glColor3f(0.5f, 0.5f, 1.0f); break;
            default:
                glColor3f(0.5f, 0.5f, 0.5f); break;               
        }
        //GLfloat gf1 = 0;
        ////GLfloat gf2, gf3 = (abc->size()%5)/5;
        //float a = 0.5+ (double(abc->size() % 5)) / 5;
        //std::cout << "a" << a << std::endl;
        //GLfloat gf2 = a;
        //GLfloat gf3 = a;
        ////std::cout << "abc->size" << abc->size() << std::endl;
        //std::cout << "gf1 gf2 gf3" << gf1 <<" "<< gf2 <<" "<< gf3 << std::endl;
        //glColor3f(gf1, gf2, gf3);
        for (i = 0; i < abc->size(); i++) {            
            glVertex3dv((*abc)[i]->Start()->Position().data());
            glVertex3dv((*abc)[i]->End()->Position().data());
        }
    }
    glEnd();
}
void DrawSkeletons()
{
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    //glColor3f(0.4f, 0.4f, 1.0f);
    glColor4f(1.0f, 0.0f, 0.f,0.5);
    glPointSize(1.5f);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < origianlPosition.size(); i++) {
        glVertex3dv(origianlPosition[i].data());
    }
    glEnd();

    glColor4f(0.0f, 0.0f, 1.f, 0.5);
    glPointSize(2.5f);
    glBegin(GL_POINTS);
    for (size_t i = 0; i < mesh.vList.size(); i++) {
        glVertex3dv(mesh.vList[i]->Position().data());
    }
    glEnd();
    glDisable(GL_LIGHTING);
}
// Wireframe render function
void DrawWireframe() {
    HEdgeList heList = mesh.Edges();
    HEdgeList bheList = mesh.BoundaryEdges();
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINES);
    size_t i;
    for (i = 0; i < heList.size(); i++) {
        glVertex3dv(heList[i]->Start()->Position().data());
        glVertex3dv(heList[i]->End()->Position().data());
    }
    glColor3f(1.0f, 0.0f, 0.0f);
    for (i = 0; i < bheList.size(); i++) {
        glVertex3dv(bheList[i]->Start()->Position().data());
        glVertex3dv(bheList[i]->End()->Position().data());
    }
    glEnd();
}

// Hidden Line render function
void DrawHiddenLine() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_FLAT);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glColor3f(0, 0, 0);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < fList.size(); i++) {
        Face* f = fList[i];
        const Eigen::Vector3d& pos1 = f->HalfEdge()->Start()->Position();
        const Eigen::Vector3d& pos2 = f->HalfEdge()->End()->Position();
        const Eigen::Vector3d& pos3 = f->HalfEdge()->Next()->End()->Position();
        glVertex3dv(pos1.data());
        glVertex3dv(pos2.data());
        glVertex3dv(pos3.data());
    }
    glEnd();
    glDisable(GL_POLYGON_OFFSET_FILL);

    DrawWireframe();
}

// Flat Shaded render function
void DrawFlatShaded() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_FLAT);
    glEnable(GL_LIGHTING);
    glColor3f(0.4f, 0.4f, 1.0f);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < fList.size(); i++) {
        Face* f = fList[i];
        const Eigen::Vector3d& pos1 = f->HalfEdge()->Start()->Position();
        const Eigen::Vector3d& pos2 = f->HalfEdge()->End()->Position();
        const Eigen::Vector3d& pos3 = f->HalfEdge()->Next()->End()->Position();
        Eigen::Vector3d normal = (pos2 - pos1).cross(pos3 - pos1);
        normal.normalize();
        glNormal3dv(normal.data());
        glVertex3dv(pos1.data());
        glVertex3dv(pos2.data());
        glVertex3dv(pos3.data());
    }
    glEnd();
    glDisable(GL_LIGHTING);
}

// Smooth Shaded render function
void DrawSmoothShaded() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    //glColor3f(0.4f, 0.4f, 1.0f);
    glColor3f(0.4f, 0.8f, 0.4f);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < fList.size(); i++) {
        Face* f = fList[i];
        Vertex* v1 = f->HalfEdge()->Start();
        Vertex* v2 = f->HalfEdge()->End();
        Vertex* v3 = f->HalfEdge()->Next()->End();

        glNormal3dv(v1->Normal().data());
        glVertex3dv(v1->Position().data());
        glNormal3dv(v2->Normal().data());
        glVertex3dv(v2->Position().data());
        glNormal3dv(v3->Normal().data());
        glVertex3dv(v3->Position().data());
    }
    glEnd();
    glDisable(GL_LIGHTING);
}

void DrawColorSmoothShaded() {
    FaceList fList = mesh.Faces();
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glColor3f(0.4f, 0.4f, 1.0f);
    glBegin(GL_TRIANGLES);
    for (size_t i = 0; i < fList.size(); i++) {
        Face* f = fList[i];
        Vertex* v1 = f->HalfEdge()->Start();
        Vertex* v2 = f->HalfEdge()->End();
        Vertex* v3 = f->HalfEdge()->Next()->End();
        glNormal3dv(v1->Normal().data());
        glColor3dv(v1->Color().data());
        glVertex3dv(v1->Position().data());
        glNormal3dv(v2->Normal().data());
        glColor3dv(v2->Color().data());
        glVertex3dv(v2->Position().data());
        glNormal3dv(v3->Normal().data());
        glColor3dv(v3->Color().data());
        glVertex3dv(v3->Position().data());
    }
    glEnd();
    glDisable(GL_LIGHTING);
}
void DrawParameterResult() {
	HEdgeList heList = mesh.Edges();
	HEdgeList bheList = mesh.BoundaryEdges();
	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_LINES);
	size_t i;
	for (i = 0; i < heList.size(); i++) {
		glVertex2d(heList[i]->Start()->GetParameterPosition()[0], heList[i]->Start()->GetParameterPosition()[1]);
		glVertex2d(heList[i]->End()->GetParameterPosition()[0], heList[i]->End()->GetParameterPosition()[1]);
	}
	glColor3f(1.0f, 0.0f, 0.0f);
	for (i = 0; i < bheList.size(); i++) {
		glVertex2d(bheList[i]->Start()->GetParameterPosition()[0], bheList[i]->Start()->GetParameterPosition()[1]);
		glVertex2d(bheList[i]->End()->GetParameterPosition()[0], bheList[i]->End()->GetParameterPosition()[1]);
	}
	glEnd();
	return;
	if (m_texture.empty())
	{
		std::cout << "not ready" << std::endl;
		return;
	}

	else
	{
		glShadeModel(GL_SMOOTH);
		glEnable(GL_LIGHTING);
		FaceList fList = mesh.Faces();
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, m_textureID);
		glBegin(GL_TRIANGLES);
		glColor3f(1.0f, 1.0f, 1.0f);
		for (size_t i = 0; i < fList.size(); i++) {
			Face* f = fList[i];
			Vertex* v1 = f->HalfEdge()->Start();
			Vertex* v2 = f->HalfEdge()->End();
			Vertex* v3 = f->HalfEdge()->Next()->End();
			glTexCoord2fv(v1->GetParameterPosition());
			//glNormal3dv(v1->Normal().data());
			glVertex2d(v1->GetParameterPosition()[0], v1->GetParameterPosition()[1]);
			glTexCoord2fv(v2->GetParameterPosition());
			//glNormal3dv(v2->Normal().data());
			glVertex2d(v2->GetParameterPosition()[0], v2->GetParameterPosition()[1]);
			glTexCoord2fv(v3->GetParameterPosition());
			//glNormal3dv(v3->Normal().data());
			glVertex2d(v3->GetParameterPosition()[0], v3->GetParameterPosition()[1]);
		}
		glEnd();
		//glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
	}


}
void DrawTexturesShaded()
{
	
	if (m_texture.empty())
	{
		std::cout << "not ready" << std::endl;
		return;
	}
	
	else
	{
		glShadeModel(GL_SMOOTH);
		glEnable(GL_LIGHTING);
		FaceList fList = mesh.Faces();
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, m_textureID);
		glBegin(GL_TRIANGLES);
		glColor3f(1.0f, 1.0f, 1.0f);
		for (size_t i = 0; i < fList.size(); i++) {
			Face* f = fList[i];
			Vertex* v1 = f->HalfEdge()->Start();
			Vertex* v2 = f->HalfEdge()->End();
			Vertex* v3 = f->HalfEdge()->Next()->End();
			glTexCoord2fv(v1->GetParameterPosition());
			//glNormal3dv(v1->Normal().data());
			glVertex3dv(v1->Position().data());
			glTexCoord2fv(v2->GetParameterPosition());
			//glNormal3dv(v2->Normal().data());
			glVertex3dv(v2->Position().data());
			glTexCoord2fv(v3->GetParameterPosition());
			//glNormal3dv(v3->Normal().data());
			glVertex3dv(v3->Position().data());
		}
		glEnd();
		//glDisable(GL_TEXTURE_2D);
		glDisable(GL_LIGHTING);
	}


	
}
void DrawSelectionRect() {
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0, winWidth, 0, winHeight);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glColor3f(1.0f, 1.0f, 1.0f);
    glRectd(downX, winHeight - downY, lastX, winHeight - lastY);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}


// draw the selected ROI vertices on the mesh
void DrawSelectedVertices() {
    VertexList vList = mesh.Vertices();
    glColor3f(1.0, 0.0, 0.0);
    glPointSize(2.0);
    glBegin(GL_POINTS);
    size_t i;
    for (i = 0; i < vList.size(); i++) {
        if (vList[i]->Flag()) {
            switch (vList[i]->Flag() % 3) {
            case 0: // handle vertices
                glColor3f(1.0, 0.3, 0.3);
                break;
            case 1:
                glColor3f(0.3, 1.0, 0.3);
                break;
            case 2:
                glColor3f(0.3, 0.3, 1.0);
                break;
            }
            glVertex3dv(vList[i]->Position().data());
        }
    }
    glEnd();
}


// GLUT keyboard callback function
void KeyboardFunc(unsigned char ch, int x, int y) {


	int IsUniform; 
	int IterNum;
    char IsK;
	char* filename=nullptr;
	switch (ch) {
	case 'o':
	case 'O':			
		m_texture = cv::imread(GetFileName());
		cv::flip(m_texture, m_texture, 0);

		cv::namedWindow("ͼƬ");
		cv::imshow("ͼƬ", m_texture);
		glDeleteTextures(1, &m_textureID);
		m_textureID = 0;

		glEnable(GL_TEXTURE_2D);
		glGenTextures(1, &m_textureID);
		glBindTexture(GL_TEXTURE_2D, m_textureID);


			// set the texture wrapping parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// set texture wrapping to GL_REPEAT (default wrapping method)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		// set texture filtering parameters
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		// image may align
		glPixelStorei(GL_UNPACK_ALIGNMENT, (m_texture.step & 3) ? 1 : 4);
		glPixelStorei(GL_UNPACK_ROW_LENGTH, m_texture.step / m_texture.elemSize());
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_texture.cols, m_texture.rows, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, m_texture.ptr());
		gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, m_texture.cols, m_texture.rows, GL_BGR_EXT, GL_UNSIGNED_BYTE, m_texture.ptr());
	
		glDisable(GL_TEXTURE_2D);

		if (m_parameterization != nullptr) delete m_parameterization;
		m_parameterization = new Parameterization(mesh);
		m_parameterization->calculate();
		break;
    case 'u':
    case 'U':
        /************************************************************************/
        /* activate the following code if you finish the corresponding functions*/
        
		std::cout << "Choose mode for explicit Umbrellasmooth (1 for Uniform, 2 for cotangent): ";
		std::cin >> IsUniform;
		std::cout << "input the number of iteration: ";
		std::cin >> IterNum;
		if (IsUniform == 1){
			for (int i = 0; i < IterNum; i++)
				mesh.UmbrellaSmooth(true);
		}
		else if (IsUniform == 2){
			for (int i = 0; i < IterNum; i++)
				mesh.UmbrellaSmooth(false);
		}

        mesh.ComputeVertexNormals();
        mesh.ComputeVertexCurvatures();
        /************************************************************************/
        break;

    case 's':
    case 'S':
		std::cout << "Choose mode for implicit Umbrellasmooth (1 for Uniform, 2 for contangent):";
		std::cin >> IsUniform;
		std::cout << "input the number of iteration: ";
		std::cin >> IterNum;
		if (IsUniform == 1){
			for (int i = 0; i < IterNum; i++)
				mesh.ImplicitUmbrellaSmooth(true);
		}
		else if (IsUniform == 2){
			for (int i = 0; i < IterNum; i++)
				mesh.ImplicitUmbrellaSmooth(false);
		}
        mesh.ComputeVertexNormals();
        mesh.ComputeVertexCurvatures();
        break;
    case 'k':
    case 'K': 
        std::cout << " *****Skeleton extraction model***** \n";
        //mesh.SkeletonParamInit();
        mesh.Skeleton();
        mesh.ComputeVertexNormals();
        mesh.ComputeVertexCurvatures();
        break;
        //std::cout << " *****Skeleton extraction model***** \n" << "Press 'K' to start:";
        //std::cin >>IsK;
        //mesh.SkeletonParamInit(); //compute the original one ring area 
        //if (IsK == 'K' || IsK =='k') {
        //    bool IsContractionFinished =false;
        //    int i = 0;
        //    while (!IsContractionFinished && (IsK == 'K' || IsK == 'k')) {
        //        i++;
        //        IsContractionFinished = mesh.Skeleton(); 
        //        std::cout << "Mesh contraction step:" << i << "Press K to continue (ESC to eixt):";
        //        std::cin >> IsK;
        //    }
        //    std::cout << "mesh contraction finished" << std::endl;
        //}
        //else
        //    std::cout << "A Wrong KeyBoard Value, Press 'K' try again" << std::endl;
        //break;


    case '1': // key '1'
        currentMode = Viewing;
        std::cout << "Viewing mode" << std::endl;
        break;
    case '2': // key '2'
        currentMode = Selection;
        std::cout << "Selection mode" << std::endl;
        break;
    case '3': // key '3'
        currentMode = Moving;
        std::cout << "Moving mode" << std::endl;
        break;
    case 27: // Esc char
        exit(0);
        break;
    }
    glutPostRedisplay();
}


// GLUT mouse callback function
void MouseFunc(int button, int state, int x, int y) {
    if (state == GLUT_DOWN) {
        downX = x;
        downY = y;
    }
    lastX = x;
    lastY = y;
    leftDown = (button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN);
    middleDown = (button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN);
    shiftDown = (glutGetModifiers() & GLUT_ACTIVE_SHIFT);

    if (currentMode == Selection && state == GLUT_UP) {
        SelectVertexByRect();
        downX = downY = lastX = lastY = 0;
        glutPostRedisplay();
        mesh.GroupingVertexFlags();
    }

    if (currentMode == Moving && state == GLUT_DOWN) {
        selectedHandleIndex = StartMove();
        std::cout << "handle " << selectedHandleIndex << " selected\n";
    }

    if (currentMode == Moving && state == GLUT_UP) {
        if (deformer) {
            // use the deformer to get a natural deformation
            deformer->Deform();
            mesh.ComputeVertexNormals();
            mesh.ComputeVertexCurvatures();
            glutPostRedisplay();
        }
    }

}


// GLUT mouse motion callback function
void MotionFunc(int x, int y) {
    switch (currentMode) {
    case Viewing:
        if (leftDown) {
            if (!shiftDown) { // rotate
                sphi += (double)(x - lastX) / 4.0;
                stheta += (double)(lastY - y) / 4.0;
            } else { // pan with shift key
                xpan += (double)(x - lastX) * sdepth / zNear / winWidth;
                ypan += (double)(lastY - y) * sdepth / zNear / winHeight;
            }
        }
        // scale
        if (middleDown)
            sdepth += (double)(lastY - y) / 10.0;
        break;

    case Selection:
        break;

    case Moving:
        GLProjector projector;
        VertexList vList = mesh.Vertices();
        int diffX = x - lastX;
        int diffY = y - lastY;
        if (diffX == 0 && diffY == 0)
            break;
        Eigen::Vector3d offset1 = projector.UnProject(x, -y, 0);
        Eigen::Vector3d offset2 = projector.UnProject(lastX, -lastY, 0);
        Eigen::Vector3d offset = (offset1 - offset2) * 20;
        for (size_t i = 0; i < vList.size(); i++) {
            if (vList[i]->Flag() == selectedHandleIndex)
                vList[i]->SetPosition(vList[i]->Position() + offset);
        }
		if (deformer) {
		}
        break;
    }

    lastX = x;
    lastY = y;
    glutPostRedisplay();
}

// select the ROI
void SelectVertexByRect() {

	/*************************/
	/* insert your code here */
	/*************************/

    // get the selection rectangle
    int x1 = downX, x2 = lastX, y1 = winHeight - downY, y2 = winHeight - lastY;
    if (x2 < x1) {
        int tmp = x1;
        x1 = x2;
        x2 = tmp;
    }
    if (y2 < y1) {
        int tmp = y1;
        y1 = y2;
        y2 = tmp;
    }
    GLProjector projector;
   
}


// perform a single editing operation
int StartMove() {
    GLProjector projector;
    VertexList vList = mesh.Vertices();

    double minDis = 0;
    int minIndex = -1;

    // find the closest vertex to the mouse position ?
    for (size_t i = 0; i < vList.size(); i++) {
        if (vList[i]->Flag())
            continue;

        Eigen::Vector3d v = projector.Project(vList[i]->Position());
        double diffX = v(0) - lastX;
        double diffY = v(1) - (winHeight - lastY);
        double dis = diffX * diffX + diffY * diffY;
        if (dis < minDis || minIndex == -1) {
            minDis = dis;
            minIndex = (int)i;
        }
    }

    return (minIndex != -1) ? vList[minIndex]->Flag() : 0;
}


// main function
void main(int argc, char** argv) {
    glutInit(&argc, argv);
    InitGL();
    InitMenu();
    if (argc >= 2)
        mesh.LoadObjFile(argv[1]);
    else
        InitGeometry();
    for (int i = 0; i < mesh.vList.size(); i++) {
        origianlPosition.push_back(mesh.vList[i]->Position());
    }
    SetBoundaryBox(mesh.MinCoord(), mesh.MaxCoord());
    mesh.SkeletonParamInit();
    mesh.ComputeVertexNormals();
    mesh.ComputeVertexCurvatures();

    /************************************************************************/
    /* activate the following code if you finish the corresponding functions*/
    mesh.DisplayMeshInfo();
    /************************************************************************/
    //DrawBoundaryLoops();
    glutMainLoop();
}


