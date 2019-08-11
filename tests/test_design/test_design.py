from primer3plus.design import Design


def test_init():
    design = Design()


def test_set(gfp):
    design = Design()
    design.set.template(gfp)
    design.set.left_sequence(gfp[60:85])
    design.set.as_generic_task()
    results = design.run()
    print(results)
    assert results
