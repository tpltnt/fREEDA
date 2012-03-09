TEMPLATE = app

TARGET = 

TSRC = Enter Trilinos source path here
TLIB = Enter Trilinos library path here

DEPENDPATH += . components diagrams dialogs paintings

INCLUDEPATH += . components diagrams paintings dialogs \
               $$TSRC/sacado/src \
               $$TSRC/sacado/src/mpl \
               $$TLIB/sacado/src \
               $$TSRC/epetra/src \
               $$TLIB/epetra/src \
               $$TSRC/teuchos/src \
               $$TLIB/teuchos/src \
               $$TSRC/amesos/src \
               $$TSRC/nox/src-epetra \
               $$TSRC/nox/src \
               $$TLIB/nox/src \
               ../libs/fftw-2.1.5/fftw \
               ../libs/fftw-2.1.5/rfftw \
               ../simulator/inout \
               ../simulator/elements \

QT += qt3support

LIBS += -L./components -lcomponents -L./diagrams -ldiagrams \
        -L./dialogs -ldialogs -L./paintings -lpaintings \
        -L$$TLIB/teuchos/src -L$$TLIB/sacado/src -lsacado \ 
        -L$$TLIB/epetra/src -lepetra -L$$TLIB/nox/src -lnox \
        -L$$TLIB/nox/src-epetra -lnoxepetra -L$$TLIB/amesos/src \
        -L$$TLIB/ifpack/src -L$$TLIB/ml/src -L$$TLIB/epetraext/src \
        -L$$TLIB/aztecoo/src -L../simulator/lib \
        -lanalysis -lio -lnetwork -lcompat -lnoxepetra -lifpack \
        -lamesos -lml -laztecoo -lepetra -lepetraext -lnox -lteuchos -lsacado \
        -llapack -lrfftw -lfftw -lblas -lgfortran -lm

# Input
HEADERS += config.h \
           element_gui.h \
           gui.h \
           guiactions.h \
           guidoc.h \
           guifile.h \
           guiinit.h \
           guimain.h \
           guiview.h \
           interface.h \
           node.h \
           viewpainter.h \
           wire.h \
           wirelabel.h \
           components/component.h \
           components/componentdialog.h \
           components/components.h \
           components/dot_bottom.h \
           components/dot_model.h \
           components/dot_top.h \
           components/equation.h \
           components/ground.h \
           components/localreference.h \
           components/plot.h \
           diagrams/curvediagram.h \
           diagrams/diagram.h \
           diagrams/diagramdialog.h \
           diagrams/diagrams.h \
           diagrams/graph.h \
           diagrams/marker.h \
           diagrams/markerdialog.h \
           diagrams/polardiagram.h \
           diagrams/psdiagram.h \
           diagrams/rect3ddiagram.h \
           diagrams/rectdiagram.h \
           diagrams/smithdiagram.h \
           diagrams/tabdiagram.h \
           dialogs/guisettingsdialog.h \
           dialogs/labeldialog.h \
           dialogs/messagebox.h \
           dialogs/newprojdialog.h \
           dialogs/settingsdialog.h \
           dialogs/simmessage.h \
           paintings/arrow.h \
           paintings/arrowdialog.h \
           paintings/ellipse.h \
           paintings/ellipsearc.h \
           paintings/filldialog.h \
           paintings/graphicline.h \
           paintings/graphictext.h \
           paintings/graphictextdialog.h \
           paintings/id_dialog.h \
           paintings/id_text.h \
           paintings/painting.h \
           paintings/paintings.h \
           paintings/portsymbol.h \
           paintings/rectangle.h 
SOURCES += element_gui.cpp \
           fREEDA_analysis.cpp \
           fREEDA_logic.cpp \
           fREEDA_lumped.cpp \
           fREEDA_nonlinear.cpp \
           fREEDA_other.cpp \
           fREEDA_sources.cpp \
           fREEDA_thermal.cpp \
           gui.cpp \
           guiactions.cpp \
           guidoc.cpp \
           guifile.cpp \
           guiinit.cpp \
           guimain.cpp \
           guiview.cpp \
           interface.cpp \
           node.cpp \
           viewpainter.cpp \
           wire.cpp \
           wirelabel.cpp 
