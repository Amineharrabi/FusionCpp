#version 330 core
uniform vec4 particleColor;
out vec4 FragColor;

void main()
{
    FragColor = particleColor;
}
