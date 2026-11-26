"""SVG primitives and CPU rasterization of the EXACT final SVG, via QtSvg."""
import os
from pathlib import Path
from xml.etree import ElementTree as ET

NS="http://www.w3.org/2000/svg"
ET.register_namespace("",NS)
_app=None
_fonts_loaded=False


def element(tag,**attrs):
    return ET.Element(f"{{{NS}}}{tag}",{k.replace("_","-"):str(v) for k,v in attrs.items()})


def label(parent,text,x,y,size=14,**attrs):
    node=element("text",x=f"{x:.2f}",y=f"{y:.2f}",font_size=size,font_family="Arial, sans-serif",**attrs)
    node.text=str(text);parent.append(node);return node


def qt_app():
    global _app,_fonts_loaded
    os.environ.setdefault("QT_QPA_PLATFORM","offscreen")
    from PySide6.QtGui import QGuiApplication
    _app=QGuiApplication.instance() or QGuiApplication([])
    if not _fonts_loaded:
        from PySide6.QtGui import QFont,QFontDatabase
        # The offscreen Windows Qt plugin may expose zero system fonts. Load
        # installed fonts explicitly; never ship copied commercial font files.
        for path in [Path('C:/Windows/Fonts/arial.ttf'),Path('C:/Windows/Fonts/msyh.ttc'),
                     Path('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')]:
            if path.is_file():QFontDatabase.addApplicationFont(str(path))
        families=QFontDatabase.families()
        if not families:raise ValueError('No usable local fonts; refusing tofu chemical labels')
        family='Arial' if 'Arial' in families else families[0]
        QFont.insertSubstitution('sans-serif',family);QFont.insertSubstitution('Arial, sans-serif',family)
        _app.setFont(QFont(family));_fonts_loaded=True
    return _app


def text_box(text,x,y,font_size):
    qt_app()
    from PySide6.QtGui import QFont,QFontMetricsF
    font=QFont("Arial");font.setPixelSize(max(1,round(font_size)))
    metrics=QFontMetricsF(font)
    return [x-1,y-metrics.ascent()-1,x+metrics.horizontalAdvance(text)+2,y+metrics.descent()+1]


def write_svg(root,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    ET.ElementTree(root).write(path,encoding="utf-8",xml_declaration=True)


def svg_to_png(svg_path:Path,png_path:Path):
    qt_app()
    from PySide6.QtSvg import QSvgRenderer
    from PySide6.QtGui import QImage,QPainter
    renderer=QSvgRenderer(str(svg_path.resolve()))
    if not renderer.isValid():raise ValueError("QtSvg rejected final SVG")
    size=renderer.defaultSize()
    if size.width()*size.height()>16000000:raise ValueError("Raster image exceeds public pixel budget")
    image=QImage(size,QImage.Format_ARGB32);image.fill(0xffffffff)
    painter=QPainter(image);renderer.render(painter);painter.end()
    if not image.save(str(png_path)):raise ValueError("Final SVG rasterization failed")
    return {"rasterizer":"QtSvg","input":"final SVG bytes","width":size.width(),"height":size.height()}
