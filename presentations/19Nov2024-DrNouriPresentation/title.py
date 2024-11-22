from manim import *

class TitleWithImage(Scene):
    def construct(self):
        
        config.frame_width = 16  # Set new width
        config.frame_height = 9  # Set new height
        self.camera.frame_width = 16
        self.camera.frame_height = 9

        height = config.frame_height
        width = config.frame_width
        # Background Image
        bg_image = ImageMobject("../images/adrien-converse-kCrrUx7US04-unsplash.jpg")  # Path to your image
        bg_image.scale_to_fit_width(config.frame_width).move_to(ORIGIN)
        self.add(bg_image)
        # Title and Subtitle Text Box
        
        title_box = Rectangle(
            height=height/2, width=width/2, color=BLACK, fill_opacity=1
        ).move_to(height*DOWN/2 + width * RIGHT/2, aligned_edge=RIGHT+DOWN)

        # verLine = Line(8* LEFT,8 *RIGHT)
        # horLine = Line(4.5 * UP,4.5* DOWN)
        # self.add(verLine,horLine)
        
        title = Text("3+1 \nDecomposition", font_size=512, color=WHITE)
        subTitle = Text("An Introduction to \nHamiltonian Formulation \nof General Relativity.", font_size=256, line_spacing=0.75, color=WHITE)
        
        title.move_to(DOWN * 0.5 + 0.5 *RIGHT, aligned_edge=LEFT+UP)
        subTitle.move_to(DOWN * 2.2 + RIGHT * 0.5, aligned_edge=LEFT+UP)
        self.add(title_box)
        self.play(Write(title),run_time=5)
        self.wait(0.5)
        self.play(Write(subTitle), run_time=5)


