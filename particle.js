class Particle{
    constructor(position){
        this.position = position;
        this.prevPosition = position.Cpy();
        this.velocity = Vector2.Zero();
        this.color = "#28b0ff";
    }
}