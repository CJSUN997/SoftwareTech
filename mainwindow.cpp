#include "mainwindow.h"
#include "ui_mainwindow.h"
#include<QFileDialog>
#include<QDebug>
#include<QMessageBox>
#include <opencv2/opencv.hpp>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->statusBar()->addPermanentWidget( ui->labellogo,1);


    ui->statusmsg->setText("Debuging!!");                  //设置label内容
    this->statusBar()->addWidget(ui->statusmsg);          //将label添加至状态栏左侧
//开发者：孙长江 张森 /n版本 V1.0 /n 版权所有：交通大学
    connect(ui->actionimport,&QAction::triggered,this,&MainWindow::ChangeImportPath);

    connect(ui->actionabout,&QAction::triggered,this,[=]{
        QMessageBox::information(this,"about","author: 孙长江 Zhang-Sen\n version:v1.0\n authority:SJTU");
    });

//    cv::Mat image = cv::imread("str");
    connect(ui->actionSmoothShaded,&QAction::triggered,this,[=]{
        ui->widget->changeintocolorsmoothshaded();
        ui->statusmsg->setText("display mode: colorsmoothshaded");
    });
    connect(ui->actionFlatShaded,&QAction::triggered,this,[=]{
        ui->widget->changeintoflatshaded();
        ui->statusmsg->setText("display mode: FlatShaded");
    });
    connect(ui->actionHiddenLine,&QAction::triggered,this,[=]{
        ui->widget->changeintohiddenline();
        ui->statusmsg->setText("display mode: HiddenLine");
    });
    connect(ui->actionWireframe,&QAction::triggered,this,[=]{
        ui->widget->changeintowireframe();
        ui->statusmsg->setText("display mode: Wireframe");
    });
    connect(ui->comboBoxGeneration,&QComboBox::currentTextChanged,this,[=]{
        QString index=ui->comboBoxGeneration->currentText();
        qDebug()<<index;
        ui->widget->GetIteration(index.toInt());
    });
}

MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow::ChangeImportPath()
{

    ui->widget->Changeimportpath();
    ui->statusmsg->setText("import successful!");
}

