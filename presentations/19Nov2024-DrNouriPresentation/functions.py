from manim import *
from wrapper import wrapper
def create_content_scene(title, text, font_size_title=160, font_size_text=198):
    class DynamicContentScene(Scene):
        def construct(self):
            # Set frame dimensions
            config.frame_width = 16
            config.frame_height = 9
            self.camera.frame_width = 16
            self.camera.frame_height = 9

            # Black background
            self.add(Rectangle(width=16, height=9, color=BLACK, fill_opacity=1).move_to(ORIGIN))
            
            # Title
            overview_heading = Text(title, font_size=font_size_title, font="IBM Plex Mono", color=WHITE)
            self.play(FadeIn(overview_heading), run_time=0.5)
            self.wait(1)
            self.play(overview_heading.animate.shift(LEFT), run_time=0.5)
            self.wait(1)

            # Wrap text
            wrapped_text = wrapper(text, 54)
            mob = []
            for t in wrapped_text:
                mob.append(Text(t, fill_opacity=1, font_size=font_size_text, 
                                font="IBM Plex Sans", color=WHITE).move_to(ORIGIN, aligned_edge=LEFT))
            
            # Display text
            if len(mob) == 1:
                self.play(FadeIn(mob[0]))
                self.wait(1)
            else:
                for i in range(len(mob)):
                    self.play(
                        FadeIn(mob[i]),
                        mob[i].animate.shift((len(mob) - 2 * i)/len(mob) * 1 * UP), 
                        run_time=0.75
                    )
            self.wait(5)

    return DynamicContentScene

def create_title_scene(title, subtitle, title_font_size=512, subtitle_font_size=256):
    class DynamicTitleScene(Scene):
        def construct(self):
            # Set frame dimensions
            config.frame_width = 16
            config.frame_height = 9
            self.camera.frame_width = 16
            self.camera.frame_height = 9

            height = config.frame_height
            width = config.frame_width

            # Background Image (optional - you might want to pass image path as parameter)
            bg_image = ImageMobject("../images/adrien-converse-kCrrUx7US04-unsplash.jpg")
            bg_image.scale_to_fit_width(config.frame_width).move_to(ORIGIN)
            self.add(bg_image)

            # Title and Subtitle Text Box
            title_box = Rectangle(
                height=height/2, width=width/2, color=BLACK, fill_opacity=1
            ).move_to(height*DOWN/2 + width * RIGHT/2, aligned_edge=RIGHT+DOWN)

            # Create title and subtitle
            title_text = Text(title, font_size=title_font_size, color=WHITE)
            subtitle_text = Text(subtitle, font_size=subtitle_font_size, 
                                 line_spacing=0.75, color=WHITE)
            
            title_text.move_to(DOWN * 0.5 + 0.5 * RIGHT, aligned_edge=LEFT+UP)
            subtitle_text.move_to(DOWN * 2.2 + RIGHT * 0.5, aligned_edge=LEFT+UP)

            # Add and animate
            self.add(title_box)
            self.play(Write(title_text), run_time=5)
            self.wait(1)
            self.play(Write(subtitle_text), run_time=5)

    return DynamicTitleScene

def create_latex_content_scene(title, latex_text, text_description=None, 
                                title_font_size=160, latex_font_size=32, 
                                description_font_size=198):
    class LatexContentScene(Scene):
        def construct(self):
            # Set frame dimensions
            config.frame_width = 16
            config.frame_height = 9
            self.camera.frame_width = 16
            self.camera.frame_height = 9

            height = config.frame_height
            width = config.frame_width

            # Black background
            self.add(Rectangle(width=16, height=9, color=BLACK, fill_opacity=1).move_to(ORIGIN))
                            # Title
            overview_heading = Text(title, font_size=title_font_size, 
                                        font="IBM Plex Mono", color=WHITE)
            if title != "":            

                self.play(FadeIn(overview_heading), run_time=0.5)
                self.wait(1)
                self.play(overview_heading.animate.shift(1.5* LEFT), run_time=0.5)
                self.wait(4)

            # LaTeX Equation
            latex_equation = MathTex(latex_text, font_size=latex_font_size, color=WHITE)

            # Position the LaTeX equation
            latex_equation.move_to(LEFT *4)

            # Add the LaTeX equation
            self.play(Transform(overview_heading, latex_equation), run_time=3)   
            wrapped_text = wrapper(text_description, 54)
            mob = []
            for t in wrapped_text:
                mob.append(Text(t, fill_opacity=1, font_size=description_font_size, 
                                font="IBM Plex Sans", color=WHITE)
                            .move_to(ORIGIN, aligned_edge=LEFT))
            
            # Display text
            if len(mob) == 1:
                self.play(FadeIn(mob[0]))
                self.wait(1)
            else:
                for i in range(len(mob)):
                    self.play(
                        FadeIn(mob[i]),
                        mob[i].animate.shift((len(mob) - 2 * i)/len(mob) * ((1 + len(mob))/(len(mob))) * UP), 
                        run_time=0.75
                        )
            self.wait(3)

    return LatexContentScene

def create_equation_transformation_scene(
    equations, 
    title=None, 
    font_size=32, 
    color=WHITE, 
    wait_time=2, 
    title_font_size=160
):
    """
    Creates a Manim scene class for transforming a list of equations.
    
    Parameters:
    - equations: List of LaTeX equation strings
    - title: Optional title for the scene
    - font_size: Font size for the equations (default 32)
    - color: Color of the equations (default WHITE)
    - wait_time: Time to wait between transformations (default 2 seconds)
    - title_font_size: Font size for the title (default 160)
    """
    class EquationTransformationScene(Scene):
        def construct(self):
            # Only add title if provided
            if title:
                overview_heading = Text(title, 
                                        font_size=title_font_size, 
                                        font="IBM Plex Mono", 
                                        color=WHITE)
                self.play(FadeIn(overview_heading), run_time=0.5)
                self.wait(1)
                self.play(overview_heading.animate.shift(LEFT), run_time=0.5)
                self.wait(0.5)
            
            # Check if equations list is empty
            if not equations:
                return
            
            # Create equations list
            eq = [
                MathTex(eq_text, 
                        font_size=font_size, 
                        color=color)
                .move_to(ORIGIN, aligned_edge=LEFT) 
                for eq_text in equations
            ]
            
            # Display equations
            if len(eq) == 1:
                self.play(FadeIn(eq[0]))
                self.wait(wait_time)
            else:
                for i in range(len(eq)):
                    self.play(
                        FadeIn(eq[i]),
                        eq[i].animate.shift((len(eq) - 2 * i)/len(eq) * 1 * UP), 
                        run_time=0.75
                    )
                    self.wait(wait_time)
    
    return EquationTransformationScene

# Example usage:
# scene_class = create_equation_transformation_scene(
#     equations=[r"f(x) = x^2", r"f(x) = 2x^2", r"f(x) = 2x^2 + 3"],
#     title="Equation Transformations"
# )